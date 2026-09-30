"""
EMIT L2A Surface Reflectance Converter
Reads EMIT L2A Reflectance, Mask, and Uncertainty NetCDF4 files,
applies GLT-based orthorectification (no KD-tree needed), and outputs:
  1. Reflectance - good wavelength bands only (default) or all 285 bands
  2. Mask - 8 QA layers (cloud, cirrus, water, spacecraft, dilated cloud,
           AOD550, H2O, aggregate flag)
  3. Per-pixel RMSE uncertainty - single band (optional, default on)
  4. Full uncertainty cube - all bands matching reflectance (optional)

EMIT spectral range: 381-2493 nm (~285 bands, ~7.4 nm spacing)
Spatial resolution: ~60 m
"""

import os
import re
import sys
import numpy as np
import h5py
import warnings

try:
    import rasterio
    from rasterio.crs import CRS
    from rasterio.transform import Affine
except ImportError:
    print("ERROR: rasterio is required. Install with: pip install rasterio")
    sys.exit(1)

# Import spectral tools module
try:
    from spectral_tools import (
        parse_spectral_args, apply_spectral_selection,
        resample_to_sensor, nearest_band_match,
        get_good_band_mask, WATER_VAPOR_REGIONS
    )
    SPECTRAL_AVAILABLE = True
except ImportError:
    print("WARNING: spectral_tools.py not found. Spectral resampling disabled.")
    print("Place spectral_tools.py in the same directory as this script.")
    SPECTRAL_AVAILABLE = False


# ============================================================================
# Constants
# ============================================================================

# Filename patterns
# EMIT_L2A_RFL_001_20260129T155601_2602910_004.nc
# EMIT_L2A_RFLUNCERT_001_20260129T155601_2602910_004.nc
# EMIT_L2A_MASK_001_20260129T155601_2602910_004.nc
RFL_PATTERN = re.compile(
    r'^EMIT_L2A_RFL_(\d+)_(\d{8}T\d{6})_(\d+)_(\d+)\.nc$',
    re.IGNORECASE
)
UNCERT_PATTERN = re.compile(
    r'^EMIT_L2A_RFLUNCERT_(\d+)_(\d{8}T\d{6})_(\d+)_(\d+)\.nc$',
    re.IGNORECASE
)
MASK_PATTERN = re.compile(
    r'^EMIT_L2A_MASK_(\d+)_(\d{8}T\d{6})_(\d+)_(\d+)\.nc$',
    re.IGNORECASE
)

# HDF5/NetCDF4 dataset paths
RFL_DATA_PATH = 'reflectance'
UNCERT_DATA_PATH = 'reflectance_uncertainty'
MASK_DATA_PATH = 'mask'
LAT_PATH = 'location/lat'
LON_PATH = 'location/lon'
GLT_X_PATH = 'location/glt_x'
GLT_Y_PATH = 'location/glt_y'
WAVELENGTHS_PATH = 'sensor_band_parameters/wavelengths'
FWHM_PATH = 'sensor_band_parameters/fwhm'
GOOD_WL_PATH = 'sensor_band_parameters/good_wavelengths'
MASK_NAMES_PATH = 'sensor_band_parameters/mask_bands'

# Mask band names (fallback if not readable from file)
DEFAULT_MASK_NAMES = [
    'Cloud Flag', 'Cirrus Flag', 'Water Flag', 'Spacecraft Flag',
    'Dilated Cloud Flag', 'AOD550', 'H2O (g cm-2)', 'Aggregate Flag'
]


# ============================================================================
# File discovery and matching
# ============================================================================

def discover_granules(input_dir):
    """
    Scan input directory for EMIT L2A file triplets (RFL + MASK + optional UNCERT).
    Match by timestamp + orbit + scene identifiers.
    Returns list of dicts with 'rfl', 'mask', 'uncert', 'granule_id' keys.
    """
    rfl_files = {}
    uncert_files = {}
    mask_files = {}

    for filename in os.listdir(input_dir):
        if not filename.lower().endswith('.nc'):
            continue

        rfl_match = RFL_PATTERN.match(filename)
        if rfl_match:
            ver, timestamp, orbit, scene = rfl_match.groups()
            key = f"{timestamp}_{orbit}_{scene}"
            rfl_files[key] = os.path.join(input_dir, filename)
            continue

        uncert_match = UNCERT_PATTERN.match(filename)
        if uncert_match:
            ver, timestamp, orbit, scene = uncert_match.groups()
            key = f"{timestamp}_{orbit}_{scene}"
            uncert_files[key] = os.path.join(input_dir, filename)
            continue

        mask_match = MASK_PATTERN.match(filename)
        if mask_match:
            ver, timestamp, orbit, scene = mask_match.groups()
            key = f"{timestamp}_{orbit}_{scene}"
            mask_files[key] = os.path.join(input_dir, filename)

    granules = []
    for key, rfl_path in sorted(rfl_files.items()):
        granule = {
            'rfl': rfl_path,
            'mask': mask_files.get(key),
            'uncert': uncert_files.get(key),
            'granule_id': key,
        }
        if granule['mask'] is None:
            print(f"  WARNING: No mask file for granule {key}")
        granules.append(granule)

    return granules


# ============================================================================
# GLT-based orthorectification
# ============================================================================

def read_glt_and_geotransform(rfl_path, bbox=None):
    """
    Read the GLT arrays and geotransform from the reflectance file.
    Optionally clip to a bounding box (min_lat, max_lat, min_lon, max_lon).
    Returns glt_x, glt_y, geotransform (Affine), crs, output dimensions.
    """
    with h5py.File(rfl_path, 'r') as f:
        glt_x_full = f[GLT_X_PATH][:]
        glt_y_full = f[GLT_Y_PATH][:]

        # Read geotransform from root attributes
        gt = f.attrs['geotransform']
        # GDAL geotransform: [ulx, xres, xskew, uly, yskew, yres]
        full_transform = Affine(gt[1], gt[2], gt[0], gt[4], gt[5], gt[3])

        # Read spatial reference
        spatial_ref = f.attrs.get('spatial_ref', b'')
        if isinstance(spatial_ref, bytes):
            spatial_ref = spatial_ref.decode('utf-8')

        if 'EPSG' in spatial_ref and '4326' in spatial_ref:
            crs = CRS.from_epsg(4326)
        else:
            crs = CRS.from_wkt(spatial_ref) if spatial_ref else CRS.from_epsg(4326)

    full_height, full_width = glt_x_full.shape

    if bbox is not None:
        min_lat, max_lat, min_lon, max_lon = bbox
        # Pixel coordinates from geotransform
        # col = (lon - ulx) / xres, row = (lat - uly) / yres
        ulx, xres = gt[0], gt[1]
        uly, yres = gt[3], gt[5]  # yres is negative

        col_min = max(0, int((min_lon - ulx) / xres))
        col_max = min(full_width, int(np.ceil((max_lon - ulx) / xres)))
        row_min = max(0, int((max_lat - uly) / yres))   # yres < 0, so max_lat gives smaller row
        row_max = min(full_height, int(np.ceil((min_lat - uly) / yres)))

        # Sanity check
        if col_min >= col_max or row_min >= row_max:
            raise ValueError(
                f"Bounding box does not overlap with GLT grid.\n"
                f"  GLT extent: UL=({ulx:.4f}, {uly:.4f}), "
                f"LR=({ulx + full_width*xres:.4f}, {uly + full_height*yres:.4f})\n"
                f"  Requested bbox: lat [{min_lat}, {max_lat}], lon [{min_lon}, {max_lon}]"
            )

        glt_x = glt_x_full[row_min:row_max, col_min:col_max]
        glt_y = glt_y_full[row_min:row_max, col_min:col_max]

        # Adjust transform for the window offset
        new_ulx = ulx + col_min * xres
        new_uly = uly + row_min * yres
        transform = Affine(gt[1], gt[2], new_ulx, gt[4], gt[5], new_uly)

        out_height, out_width = glt_x.shape
        print(f"  GLT grid (clipped): {out_width} x {out_height} "
              f"(from {full_width} x {full_height})")
        print(f"  Bbox clip: rows [{row_min}:{row_max}], cols [{col_min}:{col_max}]")
        print(f"  New UL: ({new_ulx:.4f}, {new_uly:.4f})")
    else:
        glt_x = glt_x_full
        glt_y = glt_y_full
        transform = full_transform
        out_height, out_width = full_height, full_width
        print(f"  GLT grid: {out_width} x {out_height}")

    print(f"  Geotransform: UL=({transform.c:.4f}, {transform.f:.4f}), "
          f"res=({transform.a:.6f}, {transform.e:.6f})")

    return glt_x, glt_y, transform, crs, out_width, out_height


def apply_glt(swath_data, glt_x, glt_y):
    """
    Apply GLT lookup to orthorectify swath data onto a regular grid.

    Parameters:
        swath_data: 2D (lines, samples) or 3D (lines, samples, bands) array
        glt_x: 2D array of crosstrack indices (1-based, 0=fill)
        glt_y: 2D array of downtrack indices (1-based, 0=fill)

    Returns:
        Orthorectified 2D or 3D array with NaN for unmapped pixels
    """
    out_height, out_width = glt_x.shape
    valid = (glt_x > 0) & (glt_y > 0)

    # Convert to 0-based indices
    ix = glt_x[valid].astype(np.int64) - 1
    iy = glt_y[valid].astype(np.int64) - 1

    # Clip to swath bounds
    n_lines, n_samples = swath_data.shape[0], swath_data.shape[1]
    ix = np.clip(ix, 0, n_samples - 1)
    iy = np.clip(iy, 0, n_lines - 1)

    if swath_data.ndim == 3:
        n_bands = swath_data.shape[2]
        output = np.full((out_height, out_width, n_bands), np.nan, dtype=np.float32)
        output[valid] = swath_data[iy, ix, :]
    else:
        output = np.full((out_height, out_width), np.nan, dtype=np.float32)
        output[valid] = swath_data[iy, ix]

    return output


# ============================================================================
# Data reading
# ============================================================================

def read_reflectance(rfl_path, good_only=True, spectral_opts=None):
    """
    Read reflectance data and wavelength info.

    Parameters:
        rfl_path: path to RFL NetCDF4 file
        good_only: if True, return only good wavelength bands
        spectral_opts: dict from spectral_tools.parse_spectral_args()

    Returns:
        data: 3D array (lines, samples, bands)
        wavelengths: 1D array of center wavelengths (nm)
        fwhm: 1D array of FWHM values (nm)
        band_indices: indices of selected bands in the original 285
        resample_info: dict if sensor resampling needed, else None
    """
    resample_info = None

    with h5py.File(rfl_path, 'r') as f:
        wavelengths_all = f[WAVELENGTHS_PATH][:]
        fwhm_all = f[FWHM_PATH][:]
        good_wl = f[GOOD_WL_PATH][:]

        if not good_only:
            # Include all bands, no filtering
            band_mask = np.ones(len(wavelengths_all), dtype=bool)
            print(f"  All bands: {len(wavelengths_all)}")
        elif SPECTRAL_AVAILABLE and spectral_opts and spectral_opts['mode'] != 'default':
            # Use spectral tools for band selection
            band_mask, resample_info = apply_spectral_selection(
                wavelengths_all, good_wl, spectral_opts
            )
        else:
            # Default: good wavelengths with water vapor exclusion
            if SPECTRAL_AVAILABLE:
                band_mask = get_good_band_mask(good_wl, wavelengths_all,
                                                exclude_water=True)
            else:
                band_mask = good_wl.astype(bool)
                # Manual water vapor exclusion fallback
                for wv_min, wv_max in [(1334, 1431), (1788, 1980)]:
                    water = (wavelengths_all >= wv_min) & (wavelengths_all <= wv_max)
                    band_mask = band_mask & ~water

        band_indices = np.where(band_mask)[0]
        print(f"  Selected {len(band_indices)} of {len(wavelengths_all)} bands")

        print(f"  Reading reflectance cube...")
        data = f[RFL_DATA_PATH][:, :, band_indices]

    wavelengths = wavelengths_all[band_indices]
    fwhm = fwhm_all[band_indices]

    # Set fill/nodata to NaN
    data[data == -9999] = np.nan

    print(f"  Reflectance shape: {data.shape}")
    print(f"  Wavelength range: {wavelengths.min():.1f} - {wavelengths.max():.1f} nm")

    return data, wavelengths, fwhm, band_indices, resample_info


def read_uncertainty(uncert_path, band_indices):
    """
    Read uncertainty data for specified band indices.
    """
    with h5py.File(uncert_path, 'r') as f:
        print(f"  Reading uncertainty cube...")
        data = f[UNCERT_DATA_PATH][:, :, band_indices]

    data[data == -9999] = np.nan
    print(f"  Uncertainty shape: {data.shape}")
    return data


def read_mask(mask_path):
    """
    Read mask data and band names.
    """
    with h5py.File(mask_path, 'r') as f:
        print(f"  Reading mask layers...")
        data = f[MASK_DATA_PATH][:]

        try:
            mask_names = [n.decode('utf-8') if isinstance(n, bytes) else str(n)
                          for n in f[MASK_NAMES_PATH][:]]
        except Exception:
            mask_names = DEFAULT_MASK_NAMES

    print(f"  Mask shape: {data.shape}")
    print(f"  Mask layers: {mask_names}")
    return data, mask_names


def compute_rmse_uncertainty(uncert_data):
    """
    Compute per-pixel RMSE across all bands from the uncertainty cube.
    Returns 2D array (lines, samples).
    """
    print(f"  Computing per-pixel RMSE uncertainty...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        rmse = np.sqrt(np.nanmean(uncert_data ** 2, axis=2))
    return rmse


# ============================================================================
# ENVI header writing
# ============================================================================

def build_envi_header(envi_path, height, width, num_bands, dtype, crs, transform,
                      wavelengths=None, fwhm=None, band_names=None,
                      description='', interleave='bsq'):
    """Write a custom ENVI header with wavelength and FWHM metadata."""
    dtype_map = {
        'uint8': 1, 'int16': 2, 'int32': 3, 'float32': 4,
        'float64': 5, 'uint16': 12, 'uint32': 13,
    }

    dtype_str = np.dtype(dtype).name
    envi_dtype = dtype_map.get(dtype_str, 4)

    hdr_path = envi_path if envi_path.endswith('.hdr') else envi_path + '.hdr'

    with open(hdr_path, 'w') as f:
        f.write('ENVI\n')
        f.write(f'description = {{{description}}}\n')
        f.write(f'samples = {width}\n')
        f.write(f'lines = {height}\n')
        f.write(f'bands = {num_bands}\n')
        f.write(f'header offset = 0\n')
        f.write(f'file type = ENVI Standard\n')
        f.write(f'data type = {envi_dtype}\n')
        f.write(f'interleave = {interleave}\n')
        f.write(f'byte order = 0\n')

        if transform is not None:
            x_size = abs(transform.a)
            y_size = abs(transform.e)
            ul_x = transform.c
            ul_y = transform.f
            f.write(f'map info = {{Geographic Lat/Lon, 1, 1, {ul_x}, {ul_y}, '
                    f'{x_size}, {y_size}, WGS-84}}\n')

        if crs is not None:
            f.write(f'coordinate system string = {{{crs.to_wkt()}}}\n')

        if wavelengths is not None:
            f.write('wavelength units = Nanometers\n')
            # Write wavelengths in rows of 10 for readability
            wl_values = [f'{w:.2f}' for w in wavelengths]
            f.write('wavelength = {\n')
            for i in range(0, len(wl_values), 10):
                chunk = ', '.join(wl_values[i:i + 10])
                if i + 10 < len(wl_values):
                    f.write(f'  {chunk},\n')
                else:
                    f.write(f'  {chunk}}}\n')

        if fwhm is not None:
            fwhm_values = [f'{w:.2f}' for w in fwhm]
            f.write('fwhm = {\n')
            for i in range(0, len(fwhm_values), 10):
                chunk = ', '.join(fwhm_values[i:i + 10])
                if i + 10 < len(fwhm_values):
                    f.write(f'  {chunk},\n')
                else:
                    f.write(f'  {chunk}}}\n')

        if band_names is not None:
            bn_str = ',\n  '.join(band_names)
            f.write(f'band names = {{\n  {bn_str}}}\n')

    return hdr_path


# ============================================================================
# Output writing
# ============================================================================

def write_output(data_3d, output_path, crs, transform, band_names=None,
                 wavelengths=None, fwhm=None, description='',
                 write_geotiff=True, write_envi=True,
                      overwrite=True):
    """
    Write orthorectified data to GeoTIFF and/or ENVI.
    data_3d: array (height, width, bands) or (bands, height, width)
    Internally converts to (bands, height, width) for writing.
    """
    output_files = []

    # Ensure band-first ordering
    if data_3d.ndim == 3 and data_3d.shape[2] < data_3d.shape[0]:
        # (height, width, bands) -> (bands, height, width)
        stack = np.transpose(data_3d, (2, 0, 1)).astype(np.float32)
    elif data_3d.ndim == 2:
        stack = data_3d[np.newaxis, :, :].astype(np.float32)
    else:
        stack = data_3d.astype(np.float32)

    num_bands, height, width = stack.shape

    if write_geotiff:
        gtiff_path = output_path + '.tif'
        profile = {
            'driver': 'GTiff',
            'dtype': 'float32',
            'width': width,
            'height': height,
            'count': num_bands,
            'crs': crs,
            'transform': transform,
            'nodata': np.nan,
            'compress': 'deflate',
            'predictor': 2,
            'zlevel': 6,
            'tiled': True,
            'blockxsize': 256,
            'blockysize': 256,
        }

        print(f"    Writing GeoTIFF: {os.path.basename(gtiff_path)} "
              f"({num_bands} bands, {width}x{height})")
        if not overwrite and os.path.exists(gtiff_path):
            print(f'  [skip] {os.path.basename(gtiff_path)} already exists')
        else:
            with rasterio.open(gtiff_path, 'w', **profile) as dst:
                for i in range(num_bands):
                    dst.write(stack[i], i + 1)
                    if band_names and i < len(band_names):
                        dst.set_band_description(i + 1, band_names[i])
        output_files.append(gtiff_path)

    if write_envi:
        envi_dat = output_path + '.dat'
        envi_hdr = output_path + '.hdr'

        print(f"    Writing ENVI: {os.path.basename(envi_dat)} "
              f"({num_bands} bands, {width}x{height})")

        # Write BSQ binary
        if not overwrite and os.path.exists(envi_dat):
            print(f'  [skip] {os.path.basename(envi_dat)} already exists')
        else:
            stack.tofile(envi_dat)

            build_envi_header(
                envi_hdr,
                height=height, width=width, num_bands=num_bands,
                dtype=np.float32, crs=crs, transform=transform,
                wavelengths=wavelengths, fwhm=fwhm,
                band_names=band_names, description=description,
                interleave='bsq'
        )
        output_files.append(envi_dat)
        output_files.append(envi_hdr)

    return output_files


# ============================================================================
# Main processing
# ============================================================================

def process_granule(granule, output_dir, good_only=True,
                    include_uncertainty_cube=False, include_rmse=True,
                    spectral_opts=None, bbox=None,
                    write_geotiff=True, write_envi=True):
    """
    Process a single EMIT L2A granule.

    Parameters:
        granule: dict with 'rfl', 'mask', 'uncert', 'granule_id'
        output_dir: output directory
        good_only: use only good wavelength bands
        include_uncertainty_cube: output full per-band uncertainty
        include_rmse: output per-pixel RMSE uncertainty
        spectral_opts: dict from spectral_tools.parse_spectral_args()
        bbox: tuple (min_lat, max_lat, min_lon, max_lon) for spatial clipping
        write_geotiff, write_envi: format flags
    """
    all_outputs = []
    granule_id = granule['granule_id']
    rfl_path = granule['rfl']

    # --- Read GLT and geotransform ---
    print(f"\n  Reading GLT and geotransform...")
    glt_x, glt_y, transform, crs, out_width, out_height = \
        read_glt_and_geotransform(rfl_path, bbox=bbox)

    os.makedirs(output_dir, exist_ok=True)

    # --- Reflectance ---
    print(f"\n  === Processing Reflectance ===")
    rfl_data, wavelengths, fwhm, band_indices, resample_info = \
        read_reflectance(rfl_path, good_only, spectral_opts)

    # Handle sensor resampling before orthorectification (saves memory)
    if resample_info and SPECTRAL_AVAILABLE:
        sensor = resample_info['sensor']
        if resample_info.get('nearest', False):
            # Nearest band selection
            src_indices, band_names, band_wavelengths = \
                nearest_band_match(wavelengths, sensor)
            rfl_data = rfl_data[:, :, src_indices]
            wavelengths = wavelengths[src_indices]
            fwhm = fwhm[src_indices]
        else:
            # Gaussian convolution (or RSR if available)
            rfl_data, band_names, band_wavelengths = resample_to_sensor(
                rfl_data, wavelengths, sensor,
                rsr_dir=resample_info.get('rsr_dir')
            )
            wavelengths = band_wavelengths
            fwhm = None  # FWHM not meaningful after resampling
            band_names_override = band_names
    else:
        band_names_override = None

    print(f"  Applying GLT orthorectification...")
    rfl_ortho = apply_glt(rfl_data, glt_x, glt_y)
    del rfl_data

    # Generate band names
    if band_names_override:
        rfl_band_names = band_names_override
    else:
        rfl_band_names = [f'{w:.1f} nm' for w in wavelengths]

    # Determine output suffix for sensor-matched data
    if resample_info:
        suffix = f"_{resample_info['sensor'].upper()}"
    else:
        suffix = ""

    rfl_base = os.path.join(output_dir, f"EMIT_L2A_RFL_{granule_id}{suffix}")
    outputs = write_output(
        rfl_ortho, rfl_base, crs, transform,
        band_names=rfl_band_names,
        wavelengths=wavelengths,
        fwhm=fwhm,
        description=f'EMIT L2A Surface Reflectance{suffix} ({granule_id})',
        write_geotiff=write_geotiff, write_envi=write_envi,
            overwrite=overwrite
    )
    all_outputs.extend(outputs)
    del rfl_ortho

    # --- Mask ---
    if granule['mask']:
        print(f"\n  === Processing Mask ===")
        mask_data, mask_names = read_mask(granule['mask'])

        print(f"  Applying GLT orthorectification...")
        mask_ortho = apply_glt(mask_data, glt_x, glt_y)
        del mask_data

        mask_base = os.path.join(output_dir, f"EMIT_L2A_MASK_{granule_id}")
        outputs = write_output(
            mask_ortho, mask_base, crs, transform,
            band_names=mask_names,
            description=f'EMIT L2A Mask ({granule_id})',
            write_geotiff=write_geotiff, write_envi=write_envi,
            overwrite=overwrite
        )
        all_outputs.extend(outputs)
        del mask_ortho

    # --- Uncertainty ---
    if granule['uncert']:
        if include_rmse or include_uncertainty_cube:
            print(f"\n  === Processing Uncertainty ===")
            uncert_data = read_uncertainty(granule['uncert'], band_indices)

            # Per-pixel RMSE
            if include_rmse:
                rmse = compute_rmse_uncertainty(uncert_data)
                print(f"  Applying GLT orthorectification (RMSE)...")
                rmse_ortho = apply_glt(rmse, glt_x, glt_y)
                del rmse

                rmse_base = os.path.join(output_dir,
                                         f"EMIT_L2A_RMSE_{granule_id}")
                outputs = write_output(
                    rmse_ortho, rmse_base, crs, transform,
                    band_names=['Per-pixel RMSE Uncertainty'],
                    description=f'EMIT L2A Per-pixel RMSE Uncertainty ({granule_id})',
                    write_geotiff=write_geotiff, write_envi=write_envi,
            overwrite=overwrite
                )
                all_outputs.extend(outputs)
                del rmse_ortho

            # Full uncertainty cube
            if include_uncertainty_cube:
                print(f"  Applying GLT orthorectification (full uncertainty)...")
                uncert_ortho = apply_glt(uncert_data, glt_x, glt_y)
                del uncert_data

                uncert_band_names = [f'Uncertainty {w:.1f} nm' for w in wavelengths]
                uncert_base = os.path.join(output_dir,
                                           f"EMIT_L2A_UNCERT_{granule_id}")
                outputs = write_output(
                    uncert_ortho, uncert_base, crs, transform,
                    band_names=uncert_band_names,
                    wavelengths=wavelengths,
                    fwhm=fwhm,
                    description=f'EMIT L2A Reflectance Uncertainty ({granule_id})',
                    write_geotiff=write_geotiff, write_envi=write_envi,
            overwrite=overwrite
                )
                all_outputs.extend(outputs)
                del uncert_ortho
            else:
                del uncert_data
    else:
        if include_rmse or include_uncertainty_cube:
            print(f"\n  WARNING: No uncertainty file found for this granule.")

    return all_outputs


def process_directory(input_dir, output_dir=None, good_only=True,
                      include_uncertainty_cube=False, include_rmse=True,
                      spectral_opts=None, bbox=None,
                      write_geotiff=True, write_envi=True,
                      overwrite=True):
    """Process all EMIT L2A granules found in the input directory.

    Parameters:
        bbox: tuple (min_lat, max_lat, min_lon, max_lon) for spatial clipping,
              or None for full extent.
    """
    if output_dir is None:
        output_dir = os.path.join(input_dir, 'converted')

    print(f"Scanning for EMIT L2A files in: {input_dir}")
    granules = discover_granules(input_dir)

    if not granules:
        print("No EMIT L2A granules found.")
        print("Ensure the input directory contains EMIT_L2A_RFL_*.nc files.")
        return

    print(f"Found {len(granules)} granule(s):\n")

    all_outputs = []
    for i, granule in enumerate(granules, 1):
        print(f"{'='*70}")
        print(f"Granule {i}/{len(granules)}: {granule['granule_id']}")
        print(f"  RFL:    {os.path.basename(granule['rfl'])}")
        print(f"  MASK:   {os.path.basename(granule['mask']) if granule['mask'] else 'NOT FOUND'}")
        print(f"  UNCERT: {os.path.basename(granule['uncert']) if granule['uncert'] else 'NOT FOUND'}")
        print(f"{'='*70}")

        try:
            outputs = process_granule(
                granule, output_dir,
                good_only=good_only,
                include_uncertainty_cube=include_uncertainty_cube,
                include_rmse=include_rmse,
                spectral_opts=spectral_opts,
                bbox=bbox,
                write_geotiff=write_geotiff,
                write_envi=write_envi
            )
            all_outputs.extend(outputs)
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()

        print()

    print(f"Processing complete. {len(all_outputs)} output file(s) created in: {output_dir}")


# ============================================================================
# Command-line interface
# ============================================================================

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("EMIT L2A Surface Reflectance Converter")
        print("=" * 50)
        print(f"\nUsage: python {os.path.basename(__file__)} <input_dir> [output_dir] [options]")
        print(f"\nOutput format options:")
        print(f"  --geotiff           Output GeoTIFF only")
        print(f"  --envi              Output ENVI only")
        print(f"  (default: both GeoTIFF and ENVI)")
        print(f"\nBand selection:")
        print(f"  --all-bands         Include all 285 bands (default: good wavelengths only)")
        print(f"\nSpectral resampling (requires spectral_tools.py):")
        print(f"  --region <name>     Wavelength region: vis, vnir, nir, swir, full")
        print(f"  --bands <ranges>    Custom ranges: e.g., 450-900,2000-2400")
        print(f"  --match-sensor <s>  Resample to sensor: landsat, sentinel2, aster, modis")
        print(f"  --nearest           Use nearest band instead of Gaussian convolution")
        print(f"  --every-nth <n>     Keep every Nth band")
        print(f"  --keep-water        Don't auto-exclude water vapor bands")
        print(f"  --rsr-dir <path>    Directory with RSR CSV files for convolution")
        print(f"\nSpatial subsetting:")
        print(f"  --bbox <min_lat> <max_lat> <min_lon> <max_lon>")
        print(f"                      Clip output to bounding box")
        print(f"\nUncertainty options:")
        print(f"  --no-rmse           Skip per-pixel RMSE uncertainty output")
        print(f"  --full-uncertainty  Include full per-band uncertainty cube")
        print(f"  (default: RMSE only)")
        print(f"\nExamples:")
        print(f"  python {os.path.basename(__file__)} D:\\Data\\EMIT")
        print(f"  python {os.path.basename(__file__)} D:\\Data\\EMIT --region vnir --envi")
        print(f"  python {os.path.basename(__file__)} D:\\Data\\EMIT --match-sensor landsat")
        print(f"  python {os.path.basename(__file__)} D:\\Data\\EMIT --every-nth 5")
        print(f"  python {os.path.basename(__file__)} D:\\Data\\EMIT --bands 450-900,2000-2400")
        sys.exit(0)

    args = sys.argv[1:]
    write_geotiff = True
    write_envi = True
    good_only = True
    include_rmse = True
    include_uncertainty_cube = False
    positional = []
    spectral_opts = None
    bbox = None

    # Parse spectral options first if available
    if SPECTRAL_AVAILABLE:
        spectral_opts = parse_spectral_args(args)

    # Extract bbox before positional parsing
    for i, arg in enumerate(args):
        if arg.lower() == '--bbox' and i + 4 < len(args):
            try:
                bbox = (float(args[i+1]), float(args[i+2]),
                        float(args[i+3]), float(args[i+4]))
            except ValueError:
                print("Error: --bbox requires 4 numeric values: min_lat max_lat min_lon max_lon")
                sys.exit(1)
            break

    for arg in args:
        a = arg.lower()
        if a == '--geotiff':
            write_geotiff = True
            write_envi = False
        elif a == '--envi':
            write_geotiff = False
            write_envi = True
        elif a == '--all-bands':
            good_only = False
        elif a == '--no-rmse':
            include_rmse = False
        elif a == '--full-uncertainty':
            include_uncertainty_cube = True
        elif not a.startswith('--'):
            # Skip args that are values for -- flags
            prev_idx = args.index(arg) - 1
            if prev_idx >= 0 and args[prev_idx].startswith('--'):
                continue
            positional.append(arg)

    # Re-parse positional args cleanly
    positional = []
    skip_next = False
    for i, arg in enumerate(args):
        if skip_next:
            skip_next = False
            continue
        if arg.startswith('--'):
            # Check if this flag takes a value argument
            if arg.lower() in ('--region', '--bands', '--match-sensor',
                               '--every-nth', '--rsr-dir'):
                skip_next = True
            elif arg.lower() == '--bbox':
                skip_next = False  # handled specially: skip next 4
                # Skip the 4 bbox values
                for _ in range(4):
                    i += 1
                continue
            continue
        positional.append(arg)

    if not positional:
        print("Error: No input directory specified.")
        sys.exit(1)

    input_dir = positional[0]
    output_dir = positional[1] if len(positional) > 1 else None

    if not os.path.isdir(input_dir):
        print(f"Error: Input directory does not exist: {input_dir}")
        sys.exit(1)

    band_mode = "all 285 bands" if not good_only else "good wavelengths only"
    if spectral_opts and spectral_opts['mode'] != 'default':
        band_mode = f"{spectral_opts['mode']}: {spectral_opts.get('region') or spectral_opts.get('ranges') or spectral_opts.get('sensor') or spectral_opts.get('every_nth')}"

    print(f"Output formats: {'GeoTIFF' if write_geotiff else ''}"
          f"{' + ' if write_geotiff and write_envi else ''}"
          f"{'ENVI' if write_envi else ''}")
    print(f"Band selection: {band_mode}")
    if bbox:
        print(f"Spatial clip: lat [{bbox[0]}, {bbox[1]}], lon [{bbox[2]}, {bbox[3]}]")
    print(f"RMSE uncertainty: {'yes' if include_rmse else 'no'}")
    print(f"Full uncertainty cube: {'yes' if include_uncertainty_cube else 'no'}")
    print()

    process_directory(
        input_dir, output_dir,
        good_only=good_only,
        include_uncertainty_cube=include_uncertainty_cube,
        include_rmse=include_rmse,
        spectral_opts=spectral_opts,
        bbox=bbox,
        write_geotiff=write_geotiff,
        write_envi=write_envi
    )