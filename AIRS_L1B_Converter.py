"""
AIRS L1B Radiance Converter
Reads AIRS L1B HDF4 files (HDF-EOS2 format), extracts calibrated radiances
and channel frequencies, reprojects swath data using cKDTree, and outputs
GeoTIFF and/or ENVI files.

AIRS: 2378 spectral channels, 3.75-15.4 µm (649.6-2665.2 cm-1)
Spatial resolution: ~13.5 km at nadir, 90 cross-track pixels

Output:
  1. Radiance cube (2378 bands or user-selected subset)

Channel frequencies are read from the Vdata 'nominal_freq' field and
converted to wavelength (µm) for ENVI headers.

Dependencies: pyhdf, numpy, scipy, rasterio
"""

import os
import re
import sys
import numpy as np
import warnings

try:
    from pyhdf.SD import SD, SDC
    from pyhdf.HDF import HDF, HC
    from pyhdf.VS import VS
except ImportError:
    print("ERROR: pyhdf is required. Install with: pip install pyhdf")
    sys.exit(1)

try:
    import rasterio
    from rasterio.crs import CRS
except ImportError:
    print("ERROR: rasterio is required. Install with: pip install rasterio")
    sys.exit(1)

from scipy.spatial import cKDTree


# ============================================================================
# Constants
# ============================================================================

# Filename: AIRS.2007.11.23.127.L1B.AIRS_Rad.v5.0.0.0.G07333215651.hdf
AIRS_PATTERN = re.compile(
    r'^AIRS\.(\d{4})\.(\d{2})\.(\d{2})\.(\d+)\.L1B\.AIRS_Rad\..*\.hdf$',
    re.IGNORECASE
)

# Default resolution (AIRS is ~13.5 km at nadir)
DEFAULT_RESOLUTION_KM = 13.5
DEFAULT_RESOLUTION_M = DEFAULT_RESOLUTION_KM * 1000

# KD-tree parameters (no subsampling needed - grid is already small)
K_NEIGHBORS = 4
MAX_DISTANCE_FACTOR = 3.0
N_JOBS = -1
DEG_PER_METER = 1.0 / 111320.0

# Fill value
FILL_VALUE = -9999.0


# ============================================================================
# File discovery
# ============================================================================

def discover_airs_files(input_dir):
    """Scan for AIRS L1B HDF files."""
    files = []

    def scan_dir(d):
        for f in os.listdir(d):
            match = AIRS_PATTERN.match(f)
            if match:
                year, month, day, granule = match.groups()
                gid = f"{year}{month}{day}_{granule}"
                files.append((os.path.join(d, f), gid))

    scan_dir(input_dir)
    for entry in os.listdir(input_dir):
        subdir = os.path.join(input_dir, entry)
        if os.path.isdir(subdir):
            scan_dir(subdir)

    return sorted(files, key=lambda x: x[1])


# ============================================================================
# Data reading
# ============================================================================

def read_channel_frequencies(filepath):
    """
    Read nominal channel frequencies from HDF-EOS Vdata.
    Returns frequencies in cm-1 and wavelengths in µm.
    """
    f = HDF(filepath, HC.READ)
    vs = f.vstart()

    ref = vs.find('nominal_freq')
    if not ref:
        vs.end()
        f.close()
        raise ValueError("nominal_freq Vdata not found in file")

    vd = vs.attach(ref)
    n_channels = vd.inquire()[0]
    data = vd.read(n_channels)
    vd.detach()
    vs.end()
    f.close()

    freq = np.array([d[0] for d in data], dtype=np.float64)
    wavelengths = 10000.0 / freq  # Convert cm-1 to µm

    print(f"  Channels: {len(freq)}")
    print(f"  Frequency range: {freq.min():.2f} - {freq.max():.2f} cm-1")
    print(f"  Wavelength range: {wavelengths.min():.3f} - {wavelengths.max():.3f} µm")

    return freq, wavelengths


def read_airs_data(filepath):
    """
    Read AIRS L1B radiance and geolocation data.
    """
    # Read channel frequencies first (via Vdata)
    print(f"  Reading channel frequencies...")
    freq, wavelengths = read_channel_frequencies(filepath)

    # Read science data via SD interface
    f = SD(filepath, SDC.READ)

    print(f"  Reading radiances...")
    radiance = f.select('radiances')[:].astype(np.float32)  # (track, xtrack, channel)

    print(f"  Reading geolocation...")
    lat = f.select('Latitude')[:].astype(np.float64)
    lon = f.select('Longitude')[:].astype(np.float64)

    # Read optional datasets
    try:
        state = f.select('state')[:].astype(np.int32)
    except Exception:
        state = None

    try:
        land_frac = f.select('landFrac')[:].astype(np.float32)
    except Exception:
        land_frac = None

    attrs = f.attributes()
    f.end()

    n_track, n_xtrack, n_chan = radiance.shape
    print(f"  Dimensions: {n_track} along-track x {n_xtrack} cross-track x {n_chan} channels")
    print(f"  Lat range: {lat.min():.4f} to {lat.max():.4f}")
    print(f"  Lon range: {lon.min():.4f} to {lon.max():.4f}")

    # Handle fill values
    radiance[radiance == FILL_VALUE] = np.nan
    radiance[radiance < 0] = np.nan

    # Sort by wavelength (ascending) instead of frequency
    wl_sort = np.argsort(wavelengths)
    wavelengths_sorted = wavelengths[wl_sort]
    freq_sorted = freq[wl_sort]
    radiance_sorted = radiance[:, :, wl_sort]

    print(f"  Bands reordered by wavelength (ascending)")

    return {
        'radiance': radiance_sorted,
        'lat': lat,
        'lon': lon,
        'wavelengths': wavelengths_sorted,
        'frequencies': freq_sorted,
        'state': state,
        'land_frac': land_frac,
        'attributes': attrs,
    }


# ============================================================================
# cKDTree reprojection
# ============================================================================

def compute_output_grid(lat, lon, resolution_m, bbox=None):
    """Compute output grid, optionally clipped to bbox (min_lat, max_lat, min_lon, max_lon)."""
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

    transform = rasterio.transform.from_origin(lon_min, lat_max, res_x, res_y)
    crs = CRS.from_epsg(4326)
    grid_lon_2d, grid_lat_2d = np.meshgrid(grid_lon_1d, grid_lat_1d)

    out_w, out_h = len(grid_lon_1d), len(grid_lat_1d)
    print(f"  Output grid: {out_w} x {out_h} pixels (~{resolution_m/1000:.1f} km)")
    return grid_lon_2d, grid_lat_2d, transform, crs, out_w, out_h


def build_kdtree(lat, lon):
    """Build cKDTree — no subsampling needed for AIRS (small grid)."""
    valid = np.isfinite(lat) & np.isfinite(lon) & (lat != 0) & (lon != 0)
    coords = np.column_stack([lat[valid], lon[valid]])
    src_indices = np.where(valid.ravel())[0]
    tree = cKDTree(coords)
    print(f"  KD-tree: {len(coords):,} points (no subsampling)")
    return tree, src_indices


def resample_cube_bandwise(tree, src_indices, grid_lat, grid_lon,
                           data_cube, max_dist):
    """Resample 3D cube band-by-band for memory efficiency."""
    out_height, out_width = grid_lat.shape
    n_bands = data_cube.shape[2]

    query_points = np.column_stack([grid_lat.ravel(), grid_lon.ravel()])
    n_points = len(query_points)
    print(f"  Querying KD-tree ({n_points:,} points)...")

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

    # Pre-compute neighbor indices
    neighbor_swath_idx = np.zeros((K_NEIGHBORS, n_points), dtype=np.int64)
    for k in range(K_NEIGHBORS):
        neighbor_swath_idx[k] = src_indices[indices[:, k]]

    output = np.zeros((out_height, out_width, n_bands), dtype=np.float32)

    print(f"  Resampling {n_bands} bands...")
    for b in range(n_bands):
        if (b + 1) % 500 == 0 or b == 0 or b == n_bands - 1:
            print(f"    Band {b+1}/{n_bands}...")

        flat_band = data_cube[:, :, b].ravel()
        band_out = np.zeros(n_points, dtype=np.float32)
        band_w = np.zeros(n_points, dtype=np.float32)

        for k in range(K_NEIGHBORS):
            vals = flat_band[neighbor_swath_idx[k]]
            nan_mask = np.isnan(vals)
            vals_clean = np.where(nan_mask, 0, vals)
            w = np.where(nan_mask, 0, weights[:, k])
            band_out += vals_clean * w
            band_w += w

        band_w[band_w == 0] = np.nan
        band_out = band_out / band_w
        band_out[~valid_query] = np.nan
        output[:, :, b] = band_out.reshape(out_height, out_width)

    return output


# ============================================================================
# ENVI header
# ============================================================================

def build_envi_header(envi_path, height, width, num_bands, dtype, crs, transform,
                      wavelengths=None, band_names=None, description='',
                      interleave='bsq'):
    """Write ENVI header."""
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
            ul_x, ul_y = transform.c, transform.f
            x_size, y_size = abs(transform.a), abs(transform.e)
            f.write(f'map info = {{Geographic Lat/Lon, 1, 1, {ul_x}, {ul_y}, '
                    f'{x_size}, {y_size}, WGS-84}}\n')

        if crs is not None:
            f.write(f'coordinate system string = {{{crs.to_wkt()}}}\n')

        if wavelengths is not None:
            f.write('wavelength units = Micrometers\n')
            wl_vals = [f'{w:.4f}' for w in wavelengths]
            f.write('wavelength = {\n')
            for i in range(0, len(wl_vals), 10):
                chunk = ', '.join(wl_vals[i:i+10])
                if i + 10 < len(wl_vals):
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

def write_output(data_3d, output_path, crs, transform, band_names,
                 wavelengths=None, description='',
                 write_geotiff=True, write_envi=True,
                      overwrite=True):
    """Write orthorectified data to GeoTIFF and/or ENVI."""
    output_files = []

    # Data comes from resample as (height, width, bands)
    # Need to convert to (bands, height, width) for writing
    if data_3d.ndim == 3:
        # Always assume input is (height, width, bands) from our resampler
        stack = np.transpose(data_3d, (2, 0, 1)).astype(np.float32)
    elif data_3d.ndim == 2:
        stack = data_3d[np.newaxis, :, :].astype(np.float32)
    else:
        stack = data_3d.astype(np.float32)

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
        print(f"    Writing GeoTIFF: {os.path.basename(gtiff_path)} "
              f"({num_bands} bands, {width}x{height})")
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
        print(f"    Writing ENVI: {os.path.basename(envi_dat)} "
              f"({num_bands} bands, {width}x{height})")
        if not overwrite and os.path.exists(envi_dat):
            print(f'  [skip] {os.path.basename(envi_dat)} already exists')
        else:
            stack.tofile(envi_dat)
            build_envi_header(
                envi_hdr, height, width, num_bands, np.float32, crs, transform,
                wavelengths=wavelengths, band_names=band_names,
                description=description
        )
        output_files.append(envi_dat)
        output_files.append(envi_hdr)

    return output_files


# ============================================================================
# Main processing
# ============================================================================

def process_granule(filepath, granule_id, output_dir, resolution_m=None,
                    bbox=None, spectral_opts=None,
                    write_geotiff=True, write_envi=True):
    """Process a single AIRS L1B granule."""
    all_outputs = []

    if resolution_m is None:
        resolution_m = DEFAULT_RESOLUTION_M

    # Read data
    data = read_airs_data(filepath)

    # --- Apply spectral filtering if requested ---
    if spectral_opts is not None and spectral_opts.get('mode', 'default') != 'default':
        try:
            from spectral_tools import apply_spectral_selection
            wl_nm = data['wavelengths'] * 1000.0  # µm to nm
            good_wl = np.ones(len(wl_nm), dtype=np.int32)  # all channels "good"
            band_mask, _ = apply_spectral_selection(wl_nm, good_wl, spectral_opts)
            band_indices = np.where(band_mask)[0]
            n_orig = len(data['wavelengths'])
            data['radiance'] = data['radiance'][:, :, band_indices]
            data['wavelengths'] = data['wavelengths'][band_indices]
            data['frequencies'] = data['frequencies'][band_indices]
            print(f"  Spectral filter: {len(band_indices)} of {n_orig} channels selected")
        except ImportError:
            print("  WARNING: spectral_tools.py not found, skipping spectral filter")
        except Exception as e:
            print(f"  WARNING: Spectral filtering failed: {e}")

    # Compute grid and KD-tree
    print(f"\n  Computing output grid...")
    grid_lon, grid_lat, transform, crs, out_w, out_h = \
compute_output_grid(data['lat'], data['lon'], resolution_m, bbox=bbox)

    print(f"\n  Building KD-tree...")
    tree, src_indices = build_kdtree(data['lat'], data['lon'])

    center_lat = (data['lat'].min() + data['lat'].max()) / 2.0
    pixel_deg = resolution_m * DEG_PER_METER / np.cos(np.radians(center_lat))
    max_dist = pixel_deg * MAX_DISTANCE_FACTOR

    os.makedirs(output_dir, exist_ok=True)

    # Resample radiance
    print(f"\n  === Resampling Radiance ===")
    rad_ortho = resample_cube_bandwise(
        tree, src_indices, grid_lat, grid_lon, data['radiance'], max_dist
    )

    wl = data['wavelengths']
    band_names = [f'{w:.4f} um - {10000/w:.2f} cm-1' for w in wl]

    rad_base = os.path.join(output_dir, f"AIRS_L1B_Radiance_{granule_id}")
    outputs = write_output(
        rad_ortho, rad_base, crs, transform,
        band_names=band_names, wavelengths=wl,
        description=f'AIRS L1B Calibrated Radiance ({granule_id})',
        write_geotiff=write_geotiff, write_envi=write_envi,
            overwrite=overwrite
    )
    all_outputs.extend(outputs)

    return all_outputs


def process_directory(input_dir, output_dir=None, resolution_m=None,
                      bbox=None, spectral_opts=None,
                      write_geotiff=True, write_envi=True,
                      overwrite=True):
    """Process all AIRS L1B files in a directory."""
    if output_dir is None:
        output_dir = os.path.join(input_dir, 'converted')

    print(f"Scanning for AIRS L1B files in: {input_dir}")
    files = discover_airs_files(input_dir)

    if not files:
        print("No AIRS L1B files found.")
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

    print(f"Processing complete. {len(all_outputs)} file(s) in: {output_dir}")


# ============================================================================
# CLI
# ============================================================================

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("AIRS L1B Radiance Converter")
        print("=" * 50)
        print(f"\nUsage: python {os.path.basename(__file__)} <input_dir> [output_dir] [options]")
        print(f"\nOutput format:")
        print(f"  --geotiff         GeoTIFF only")
        print(f"  --envi            ENVI only")
        print(f"  (default: both)")
        print(f"\nResolution:")
        print(f"  --resolution <m>  Output pixel size in meters (default: {DEFAULT_RESOLUTION_M:.0f})")
        print(f"\nExamples:")
        print(f"  python {os.path.basename(__file__)} D:\\Data\\AIRS")
        print(f"  python {os.path.basename(__file__)} D:\\Data\\AIRS --envi")
        print(f"\nAIRS: 2378 channels, 3.75-15.4 µm, ~13.5 km resolution")
        print(f"Channels are output in wavelength order (ascending).")
        print(f"Band names include both wavelength (µm) and wavenumber (cm-1).")
        sys.exit(0)

    args = sys.argv[1:]
    write_geotiff = True
    write_envi = True
    resolution = None
    positional = []

    i = 0
    while i < len(args):
        a = args[i].lower()
        if a == '--geotiff':
            write_geotiff = True; write_envi = False
        elif a == '--envi':
            write_geotiff = False; write_envi = True
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

    process_directory(input_dir, output_dir, resolution_m=resolution,
                      write_geotiff=write_geotiff, write_envi=write_envi)