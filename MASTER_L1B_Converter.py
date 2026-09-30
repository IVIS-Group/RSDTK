"""
MASTER L1B Calibrated Radiance Converter
Reads MASTER L1B HDF4 files, reorders bands by wavelength, applies optional
VSWIR reflectance conversion, reprojects swath data using cKDTree, and outputs
multiband GeoTIFF and/or ENVI files.

Outputs:
  1. VSWIR (Bands 1-25): Radiance or TOA Reflectance (user choice)
  2. TIR (Bands 26-50): Radiance only

MASTER has 50 spectral bands spanning 0.46-12.81 µm.
Note: Bands are reordered by wavelength (band 26 at 4.07 µm is out of sequence).

Dependencies: pyhdf, numpy, scipy, rasterio
"""

import os
import re
import sys
import numpy as np
import warnings

try:
    from pyhdf.SD import SD, SDC
except ImportError:
    print("ERROR: pyhdf is required. Install with: pip install pyhdf")
    sys.exit(1)

try:
    import rasterio
    from rasterio.crs import CRS
    from rasterio.transform import Affine
except ImportError:
    print("ERROR: rasterio is required. Install with: pip install rasterio")
    sys.exit(1)

from scipy.spatial import cKDTree


# ============================================================================
# Constants
# ============================================================================

# Filename pattern: MASTERL1B_2598200_10_20250923_2129_2139_V01.hdf
L1B_PATTERN = re.compile(
    r'^MASTERL1B_(\d+)_(\d+)_(\d{8})_(\d{4})_(\d{4})_V(\d+)\.hdf$',
    re.IGNORECASE
)

# HDF4 dataset names
CALIBRATED_DATA = 'CalibratedData'
PIXEL_LAT = 'PixelLatitude'
PIXEL_LON = 'PixelLongitude'
CENTER_WL = 'Central100%ResponseWavelength'
LEFT_WL = 'Left50%ResponseWavelength'
RIGHT_WL = 'Right50%ResponseWavelength'
SOLAR_IRRADIANCE = 'SolarSpectralIrradiance'
SOLAR_ZENITH = 'SolarZenithAngle'

# VSWIR/TIR split: bands with center wavelength < 3.0 µm are VSWIR
VSWIR_THRESHOLD_UM = 3.0

# cKDTree parameters
SUBSAMPLE_FACTOR = 3
K_NEIGHBORS = 4
MAX_DISTANCE_FACTOR = 3.0
N_JOBS = -1
DEG_PER_METER = 1.0 / 111320.0


# ============================================================================
# File discovery
# ============================================================================

def discover_l1b_files(input_dir):
    """
    Scan for MASTER L1B HDF4 files.
    Looks in input_dir and one level of subdirectories.
    Returns list of (filepath, granule_id) tuples.
    """
    files = []

    def scan_dir(d):
        for f in os.listdir(d):
            match = L1B_PATTERN.match(f)
            if match:
                flight, line, date, start, end, ver = match.groups()
                gid = f"{flight}_{line}_{date}_{start}_{end}"
                files.append((os.path.join(d, f), gid))

    scan_dir(input_dir)
    for entry in os.listdir(input_dir):
        subdir = os.path.join(input_dir, entry)
        if os.path.isdir(subdir):
            scan_dir(subdir)

    return sorted(files, key=lambda x: x[1])


# ============================================================================
# HDF4 reading
# ============================================================================

def read_l1b_data(filepath):
    """
    Read all science data and metadata from MASTER L1B HDF4 file.

    Returns dict with:
        'radiance': (lines, bands, samples) -> reordered to (lines, samples, bands_sorted)
        'wavelengths': sorted center wavelengths (µm)
        'left_wl': sorted left band edges (µm)
        'right_wl': sorted right band edges (µm)
        'solar_irradiance': per-band solar irradiance
        'lat': (lines, samples)
        'lon': (lines, samples)
        'solar_zenith': (lines, samples) or None
        'sort_order': indices that sort bands by wavelength
        'original_band_nums': original 1-based band numbers in sorted order
        'attributes': file-level metadata dict
    """
    f = SD(filepath, SDC.READ)

    # Read datasets
    print(f"  Reading calibrated radiance...")
    cal_ds = f.select(CALIBRATED_DATA)
    cal_data = cal_ds[:]  # (lines, bands, samples)
    cal_attrs = cal_ds.attributes()
    scale_factors = np.array(cal_attrs.get('scale_factor', np.ones(50)), dtype=np.float64)
    fill_value = cal_attrs.get('_FillValue', -999)

    print(f"  Reading geolocation...")
    lat = f.select(PIXEL_LAT)[:].astype(np.float64)
    lon = f.select(PIXEL_LON)[:].astype(np.float64)

    center_wl = f.select(CENTER_WL)[:].astype(np.float64)
    left_wl = f.select(LEFT_WL)[:].astype(np.float64)
    right_wl = f.select(RIGHT_WL)[:].astype(np.float64)
    solar_irr = f.select(SOLAR_IRRADIANCE)[:].astype(np.float64)

    # Try to read solar zenith angle
    try:
        solar_zen = f.select(SOLAR_ZENITH)[:].astype(np.float32)
    except Exception:
        solar_zen = None
        print(f"  WARNING: SolarZenithAngle not found, reflectance conversion unavailable")

    # Read attributes
    attrs = f.attributes()

    f.end()

    n_lines, n_bands, n_samples = cal_data.shape
    print(f"  Data dimensions: {n_lines} lines x {n_samples} samples x {n_bands} bands")
    print(f"  Lat range: {lat[lat!=0].min():.4f} to {lat[lat!=0].max():.4f}")
    print(f"  Lon range: {lon[lon!=0].min():.4f} to {lon[lon!=0].max():.4f}")

    # Sort bands by wavelength
    sort_order = np.argsort(center_wl)
    original_band_nums = sort_order + 1  # 1-based

    # Check if reordering is needed
    if not np.array_equal(sort_order, np.arange(n_bands)):
        n_reordered = (sort_order != np.arange(n_bands)).sum()
        print(f"  Reordering {n_reordered} bands to wavelength order")

    # Reorder everything by wavelength
    center_wl_sorted = center_wl[sort_order]
    left_wl_sorted = left_wl[sort_order]
    right_wl_sorted = right_wl[sort_order]
    solar_irr_sorted = solar_irr[sort_order]

    # Reorder scale factors by wavelength
    scale_factors_sorted = scale_factors[sort_order]

    # Rearrange radiance: (lines, bands, samples) -> (lines, samples, bands_sorted)
    cal_sorted = cal_data[:, sort_order, :]
    cal_sorted = np.transpose(cal_sorted, (0, 2, 1))  # (lines, samples, bands)

    # Convert to float32, apply scale factors, handle fill values
    radiance = cal_sorted.astype(np.float32)

    # Mark fill values BEFORE scaling
    fill_mask = (cal_sorted == fill_value)

    # Apply per-band scale factors: Radiance = DN * scale_factor
    radiance = radiance * scale_factors_sorted[np.newaxis, np.newaxis, :]

    # Set fill pixels to NaN
    radiance[fill_mask] = np.nan

    print(f"  Scale factors applied (range: {scale_factors_sorted.min():.6f} - {scale_factors_sorted.max():.6f})")
    print(f"  Wavelength range: {center_wl_sorted[0]:.3f} - {center_wl_sorted[-1]:.3f} µm")

    return {
        'radiance': radiance,
        'wavelengths': center_wl_sorted,
        'left_wl': left_wl_sorted,
        'right_wl': right_wl_sorted,
        'solar_irradiance': solar_irr_sorted,
        'lat': lat,
        'lon': lon,
        'solar_zenith': solar_zen,
        'sort_order': sort_order,
        'original_band_nums': original_band_nums,
        'attributes': attrs,
    }


def split_vswir_tir(data):
    """
    Split the wavelength-sorted data into VSWIR and TIR groups.
    VSWIR: wavelength < 3.0 µm
    TIR: wavelength >= 3.0 µm
    """
    wl = data['wavelengths']
    vswir_mask = wl < VSWIR_THRESHOLD_UM
    tir_mask = ~vswir_mask

    vswir_idx = np.where(vswir_mask)[0]
    tir_idx = np.where(tir_mask)[0]

    print(f"  VSWIR bands: {len(vswir_idx)} ({wl[vswir_idx[0]]:.3f} - {wl[vswir_idx[-1]]:.3f} µm)")
    print(f"  TIR bands: {len(tir_idx)} ({wl[tir_idx[0]]:.3f} - {wl[tir_idx[-1]]:.3f} µm)")

    return vswir_idx, tir_idx


def compute_reflectance(radiance, solar_irradiance, solar_zenith):
    """
    Convert radiance to TOA reflectance.
    ρ = (π * L) / (E_sun * cos(θ_sun))

    Parameters:
        radiance: (lines, samples, bands) array
        solar_irradiance: (bands,) per-band solar irradiance
        solar_zenith: (lines, samples) solar zenith angle in degrees
    """
    print(f"  Converting VSWIR to TOA reflectance...")
    cos_sza = np.cos(np.radians(solar_zenith))[:, :, np.newaxis]
    cos_sza[cos_sza <= 0] = np.nan  # Avoid division by zero for night pixels

    irr = solar_irradiance[np.newaxis, np.newaxis, :]

    reflectance = (np.pi * radiance) / (irr * cos_sza)
    return reflectance.astype(np.float32)


# ============================================================================
# cKDTree reprojection
# ============================================================================

def compute_output_grid(lat, lon, resolution_m=50.0, bbox=None):
    """
    Compute output grid from swath lat/lon.
    MASTER resolution varies by altitude; default ~50m for ER-2.
    """
    valid = np.isfinite(lat) & np.isfinite(lon) & (lat != 0) & (lon != 0)
    if not valid.any():
        raise ValueError("No valid geolocation data found.")

    lat_min, lat_max = lat[valid].min(), lat[valid].max()
    lon_min, lon_max = lon[valid].min(), lon[valid].max()

    # Apply bbox clipping if provided
    if bbox is not None:
        b_lat_min, b_lat_max, b_lon_min, b_lon_max = bbox
        lat_min = max(lat_min, b_lat_min)
        lat_max = min(lat_max, b_lat_max)
        lon_min = max(lon_min, b_lon_min)
        lon_max = min(lon_max, b_lon_max)
        if lat_min >= lat_max or lon_min >= lon_max:
            raise ValueError(f"Bbox does not overlap with swath extent.")
        print(f"  Clipped to bbox: lat [{lat_min:.4f}, {lat_max:.4f}], "
              f"lon [{lon_min:.4f}, {lon_max:.4f}]")

    center_lat = (lat_min + lat_max) / 2.0
    res_y = resolution_m * DEG_PER_METER
    res_x = resolution_m * DEG_PER_METER / np.cos(np.radians(center_lat))

    grid_lon_1d = np.arange(lon_min, lon_max, res_x)
    grid_lat_1d = np.arange(lat_max, lat_min, -res_y)

    out_width = len(grid_lon_1d)
    out_height = len(grid_lat_1d)

    transform = rasterio.transform.from_origin(lon_min, lat_max, res_x, res_y)
    crs = CRS.from_epsg(4326)

    grid_lon_2d, grid_lat_2d = np.meshgrid(grid_lon_1d, grid_lat_1d)

    print(f"  Output grid: {out_width} x {out_height} pixels")
    print(f"  Resolution: {res_x:.6f} x {res_y:.6f} deg (~{resolution_m}m)")

    return grid_lon_2d, grid_lat_2d, transform, crs, out_width, out_height


def build_kdtree(lat, lon):
    """Build cKDTree from swath coordinates with subsampling."""
    valid = np.isfinite(lat) & np.isfinite(lon) & (lat != 0) & (lon != 0)

    subsample = np.zeros_like(valid, dtype=bool)
    subsample[::SUBSAMPLE_FACTOR, ::SUBSAMPLE_FACTOR] = True
    use_mask = valid & subsample

    lat_use = lat[use_mask]
    lon_use = lon[use_mask]
    src_indices = np.where(use_mask.ravel())[0]

    print(f"  KD-tree: {len(lat_use):,} points (subsample {SUBSAMPLE_FACTOR}x)")

    coords = np.column_stack([lat_use, lon_use])
    tree = cKDTree(coords)

    return tree, src_indices


def resample_cube(tree, src_indices, lat, lon, grid_lat, grid_lon,
                  data_cube, max_dist):
    """
    Resample a 3D cube (lines, samples, bands) onto regular grid.
    """
    out_height, out_width = grid_lat.shape
    n_bands = data_cube.shape[2]

    query_points = np.column_stack([grid_lat.ravel(), grid_lon.ravel()])
    print(f"  Querying KD-tree ({len(query_points):,} points)...")

    distances, indices = tree.query(query_points, k=K_NEIGHBORS, workers=N_JOBS)
    if K_NEIGHBORS == 1:
        distances = distances[:, np.newaxis]
        indices = indices[:, np.newaxis]

    valid_query = distances[:, 0] < max_dist

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        weights = 1.0 / distances
        weights[~np.isfinite(weights)] = 1e10

    w_sum = weights.sum(axis=1, keepdims=True)
    w_sum[w_sum == 0] = 1.0
    weights = weights / w_sum

    flat_data = data_cube.reshape(-1, n_bands)

    # Initialize output to ZERO (not NaN) for accumulation
    output = np.zeros((len(query_points), n_bands), dtype=np.float32)
    weight_accum = np.zeros((len(query_points), n_bands), dtype=np.float32)

    print(f"  Resampling {n_bands} bands...")
    for k in range(K_NEIGHBORS):
        swath_idx = src_indices[indices[:, k]]
        vals = flat_data[swath_idx]
        w = weights[:, k:k+1]  # (npts, 1)

        nan_mask = np.isnan(vals)
        vals_clean = np.where(nan_mask, 0, vals)
        w_expanded = np.broadcast_to(w, vals.shape).copy()
        w_expanded[nan_mask] = 0

        output += vals_clean * w_expanded
        weight_accum += w_expanded

    # Normalize by actual accumulated weights
    weight_accum[weight_accum == 0] = np.nan
    output = output / weight_accum

    # Apply distance mask
    output[~valid_query] = np.nan

    # Normalize by actual accumulated weights
    weight_accum[weight_accum == 0] = np.nan
    output = output / weight_accum

    # Apply distance mask
    output[~valid_query] = np.nan

    return output.reshape(out_height, out_width, n_bands)


# ============================================================================
# ENVI header writing
# ============================================================================

def build_envi_header(envi_path, height, width, num_bands, dtype, crs, transform,
                      wavelengths=None, fwhm=None, band_names=None,
                      description='', interleave='bsq', wavelength_units='Micrometers'):
    """Write ENVI header with wavelength metadata."""
    dtype_map = {
        'uint8': 1, 'int16': 2, 'int32': 3, 'float32': 4,
        'float64': 5, 'uint16': 12, 'uint32': 13,
    }
    envi_dtype = dtype_map.get(np.dtype(dtype).name, 4)
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
            f.write(f'wavelength units = {wavelength_units}\n')
            wl_str = ', '.join(f'{w:.4f}' for w in wavelengths)
            f.write(f'wavelength = {{\n  {wl_str}}}\n')

        if fwhm is not None:
            fwhm_str = ', '.join(f'{w:.4f}' for w in fwhm)
            f.write(f'fwhm = {{\n  {fwhm_str}}}\n')

        if band_names is not None:
            bn_str = ',\n  '.join(band_names)
            f.write(f'band names = {{\n  {bn_str}}}\n')

    return hdr_path


# ============================================================================
# Output writing
# ============================================================================

def write_output(data_3d, output_path, crs, transform, band_names,
                 wavelengths=None, fwhm=None, description='',
                 wavelength_units='Micrometers',
                 write_geotiff=True, write_envi=True,
                      overwrite=True):
    """Write orthorectified data to GeoTIFF and/or ENVI."""
    output_files = []

    # Ensure (bands, height, width)
    if data_3d.ndim == 3:
        if data_3d.shape[2] < data_3d.shape[0]:
            stack = np.transpose(data_3d, (2, 0, 1)).astype(np.float32)
        else:
            stack = data_3d.astype(np.float32)
    else:
        stack = data_3d[np.newaxis, :, :].astype(np.float32)

    num_bands, height, width = stack.shape

    if write_geotiff:
        gtiff_path = output_path + '.tif'
        profile = {
            'driver': 'GTiff', 'dtype': 'float32',
            'width': width, 'height': height, 'count': num_bands,
            'crs': crs, 'transform': transform, 'nodata': np.nan,
            'compress': 'deflate', 'predictor': 2, 'zlevel': 6,
            'tiled': True, 'blockxsize': 256, 'blockysize': 256,
        }
        print(f"    Writing GeoTIFF: {os.path.basename(gtiff_path)} ({num_bands} bands)")
        if not overwrite and os.path.exists(gtiff_path):
            print(f'  [skip] {os.path.basename(gtiff_path)} already exists')
        else:
            with rasterio.open(gtiff_path, 'w', **profile) as dst:
                for i in range(num_bands):
                    dst.write(stack[i], i + 1)
                    if i < len(band_names):
                        dst.set_band_description(i + 1, band_names[i])
        output_files.append(gtiff_path)

    if write_envi:
        envi_dat = output_path + '.dat'
        envi_hdr = output_path + '.hdr'
        print(f"    Writing ENVI: {os.path.basename(envi_dat)} ({num_bands} bands)")
        if not overwrite and os.path.exists(envi_dat):
            print(f'  [skip] {os.path.basename(envi_dat)} already exists')
        else:
            stack.tofile(envi_dat)
            build_envi_header(
                envi_hdr, height, width, num_bands, np.float32, crs, transform,
                wavelengths=wavelengths, fwhm=fwhm, band_names=band_names,
                description=description, wavelength_units=wavelength_units
        )
        output_files.append(envi_dat)
        output_files.append(envi_hdr)

    return output_files


# ============================================================================
# Main processing
# ============================================================================

def process_granule(filepath, granule_id, output_dir, vswir_as_reflectance=False,
                    resolution_m=50.0, bbox=None, spectral_opts=None,
                    
                    write_geotiff=True, write_envi=True):
    """Process a single MASTER L1B granule."""
    all_outputs = []

    # Read data
    print(f"\n  Reading L1B data...")
    data = read_l1b_data(filepath)

    # --- Apply spectral filtering if requested ---
    if spectral_opts is not None and spectral_opts.get('mode', 'default') != 'default':
        try:
            from spectral_tools import apply_spectral_selection
            wl_nm = data['wavelengths'] * 1000.0  # µm to nm
            good_wl = np.ones(len(wl_nm), dtype=np.int32)
            band_mask, _ = apply_spectral_selection(wl_nm, good_wl, spectral_opts)
            band_indices = np.where(band_mask)[0]
            n_orig = len(data['wavelengths'])
            data['radiance'] = data['radiance'][:, :, band_indices]
            data['wavelengths'] = data['wavelengths'][band_indices]
            if data.get('fwhm') is not None:
                data['fwhm'] = data['fwhm'][band_indices]
            if data.get('solar_irradiance') is not None:
                data['solar_irradiance'] = data['solar_irradiance'][band_indices]
            print(f"  Spectral filter: {len(band_indices)} of {n_orig} bands selected")
        except ImportError:
            print("  WARNING: spectral_tools.py not found, skipping spectral filter")
        except Exception as e:
            print(f"  WARNING: Spectral filtering failed: {e}")

    # Split into VSWIR and TIR
    vswir_idx, tir_idx = split_vswir_tir(data)

    # Build output grid and KD-tree
    print(f"\n  Computing output grid...")
    grid_lon, grid_lat, transform, crs, out_w, out_h = \
compute_output_grid(data['lat'], data['lon'], resolution_m, bbox=bbox)

    print(f"\n  Building KD-tree...")
    tree, src_indices = build_kdtree(data['lat'], data['lon'])

    center_lat = (data['lat'][data['lat'] != 0].min() +
                  data['lat'][data['lat'] != 0].max()) / 2.0
    pixel_deg = resolution_m * DEG_PER_METER / np.cos(np.radians(center_lat))
    max_dist = pixel_deg * MAX_DISTANCE_FACTOR

    os.makedirs(output_dir, exist_ok=True)

    # --- VSWIR ---
    print(f"\n  === Processing VSWIR ===")
    vswir_cube = data['radiance'][:, :, vswir_idx]

    if vswir_as_reflectance and data['solar_zenith'] is not None:
        vswir_cube = compute_reflectance(
            vswir_cube,
            data['solar_irradiance'][vswir_idx],
            data['solar_zenith']
        )
        vswir_label = 'TOARef'
        vswir_desc = 'TOA Reflectance'
    else:
        if vswir_as_reflectance and data['solar_zenith'] is None:
            print(f"  WARNING: Cannot compute reflectance without solar zenith, keeping radiance")
        vswir_label = 'Radiance'
        vswir_desc = 'Calibrated Radiance'

    vswir_ortho = resample_cube(
        tree, src_indices, data['lat'], data['lon'],
        grid_lat, grid_lon, vswir_cube, max_dist
    )
    del vswir_cube

    vswir_wl = data['wavelengths'][vswir_idx]
    vswir_fwhm = data['right_wl'][vswir_idx] - data['left_wl'][vswir_idx]
    vswir_bnames = [f"Band {data['original_band_nums'][i]} ({vswir_wl[j]:.3f} um)"
                    for j, i in enumerate(vswir_idx)]

    vswir_base = os.path.join(output_dir, f"MASTER_L1B_VSWIR_{vswir_label}_{granule_id}")
    outputs = write_output(
        vswir_ortho, vswir_base, crs, transform,
        band_names=vswir_bnames, wavelengths=vswir_wl, fwhm=vswir_fwhm,
        description=f'MASTER L1B VSWIR {vswir_desc} ({granule_id})',
        write_geotiff=write_geotiff, write_envi=write_envi,
            overwrite=overwrite
    )
    all_outputs.extend(outputs)
    del vswir_ortho

    # --- TIR ---
    print(f"\n  === Processing TIR ===")
    tir_cube = data['radiance'][:, :, tir_idx]

    tir_ortho = resample_cube(
        tree, src_indices, data['lat'], data['lon'],
        grid_lat, grid_lon, tir_cube, max_dist
    )
    del tir_cube

    tir_wl = data['wavelengths'][tir_idx]
    tir_fwhm = data['right_wl'][tir_idx] - data['left_wl'][tir_idx]
    tir_bnames = [f"Band {data['original_band_nums'][i]} ({tir_wl[j]:.3f} um)"
                  for j, i in enumerate(tir_idx)]

    tir_base = os.path.join(output_dir, f"MASTER_L1B_TIR_Radiance_{granule_id}")
    outputs = write_output(
        tir_ortho, tir_base, crs, transform,
        band_names=tir_bnames, wavelengths=tir_wl, fwhm=tir_fwhm,
        description=f'MASTER L1B TIR Calibrated Radiance ({granule_id})',
        write_geotiff=write_geotiff, write_envi=write_envi,
            overwrite=overwrite
    )
    all_outputs.extend(outputs)
    del tir_ortho

    return all_outputs


def process_directory(input_dir, output_dir=None, vswir_as_reflectance=False,
                      resolution_m=50.0, bbox=None, spectral_opts=None,
                      
                    write_geotiff=True, write_envi=True,
                      overwrite=True):
    """Process all MASTER L1B files in a directory."""
    if output_dir is None:
        output_dir = os.path.join(input_dir, 'converted')

    print(f"Scanning for MASTER L1B files in: {input_dir}")
    files = discover_l1b_files(input_dir)

    if not files:
        print("No MASTER L1B files found.")
        return

    print(f"Found {len(files)} file(s):\n")

    all_outputs = []
    for i, (filepath, gid) in enumerate(files, 1):
        print(f"{'='*70}")
        print(f"Granule {i}/{len(files)}: {gid}")
        print(f"  File: {os.path.basename(filepath)}")
        print(f"{'='*70}")

        try:
            outputs = process_granule(
                filepath, gid, output_dir,
                vswir_as_reflectance=vswir_as_reflectance,
                resolution_m=resolution_m,
                bbox=bbox, spectral_opts=spectral_opts,
                write_geotiff=write_geotiff, write_envi=write_envi,
            overwrite=overwrite
            )
            all_outputs.extend(outputs)
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
        print()

    print(f"Processing complete. {len(all_outputs)} file(s) created in: {output_dir}")


# ============================================================================
# CLI
# ============================================================================

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("MASTER L1B Calibrated Radiance Converter")
        print("=" * 50)
        print(f"\nUsage: python {os.path.basename(__file__)} <input_dir> [output_dir] [options]")
        print(f"\nOutput format:")
        print(f"  --geotiff         GeoTIFF only")
        print(f"  --envi            ENVI only")
        print(f"  (default: both)")
        print(f"\nVSWIR options:")
        print(f"  --reflectance     Convert VSWIR to TOA reflectance (default: radiance)")
        print(f"\nResolution:")
        print(f"  --resolution <m>  Output pixel size in meters (default: 50)")
        print(f"\nExamples:")
        print(f"  python {os.path.basename(__file__)} D:\\Data\\MASTER")
        print(f"  python {os.path.basename(__file__)} D:\\Data\\MASTER --reflectance --envi")
        print(f"  python {os.path.basename(__file__)} D:\\Data\\MASTER --resolution 25")
        print(f"\nOutputs per granule:")
        print(f"  1. VSWIR (25 bands): Radiance or TOA Reflectance")
        print(f"  2. TIR (25 bands): Calibrated Radiance")
        print(f"\nBands are automatically reordered by wavelength.")
        sys.exit(0)

    args = sys.argv[1:]
    write_geotiff = True
    write_envi = True
    vswir_reflectance = False
    resolution = 50.0
    positional = []

    i = 0
    while i < len(args):
        a = args[i].lower()
        if a == '--geotiff':
            write_geotiff = True; write_envi = False
        elif a == '--envi':
            write_geotiff = False; write_envi = True
        elif a == '--reflectance':
            vswir_reflectance = True
        elif a == '--resolution':
            resolution = float(args[i + 1]); i += 1
        else:
            positional.append(args[i])
        i += 1

    input_dir = positional[0]
    output_dir = positional[1] if len(positional) > 1 else None

    if not os.path.isdir(input_dir):
        print(f"Error: Input directory does not exist: {input_dir}")
        sys.exit(1)

    print(f"VSWIR output: {'TOA Reflectance' if vswir_reflectance else 'Radiance'}")
    print(f"Resolution: {resolution}m")
    print()

    process_directory(input_dir, output_dir,
                      vswir_as_reflectance=vswir_reflectance,
                      resolution_m=resolution,
                      write_geotiff=write_geotiff, write_envi=write_envi)