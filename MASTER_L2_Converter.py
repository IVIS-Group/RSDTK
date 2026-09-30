"""
MASTER L2 Emissivity/LST Converter
Reads MASTER L2 HDF5 files, reprojects swath data using cKDTree, and outputs:
  1. Emissivity (6 TIR bands)
  2. LST (single band, Kelvin)
  3. QA Map (single band)

Dependencies: h5py, numpy, scipy, rasterio
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
except ImportError:
    print("ERROR: rasterio is required. Install with: pip install rasterio")
    sys.exit(1)

from scipy.spatial import cKDTree


# ============================================================================
# Constants
# ============================================================================

# Filename: MASTERL2_2598200_01_20250923_1707_1730_V01_B200_SV01.hdf5
L2_PATTERN = re.compile(
    r'^MASTERL2_(\d+)_(\d+)_(\d{8})_(\d{4})_(\d{4})_V(\d+)_(\w+)_(\w+)\.hdf5$',
    re.IGNORECASE
)

# MASTER TIR emissivity band wavelengths (µm)
# These correspond to MASTER bands 42-47 (the 6 primary TIR bands used in TES)
EMIS_WAVELENGTHS_UM = [8.25, 8.65, 9.05, 9.71, 10.11, 10.58]

EMIS_BAND_NAMES = [
    'Emissivity Band 1 (8.25 um)',
    'Emissivity Band 2 (8.65 um)',
    'Emissivity Band 3 (9.05 um)',
    'Emissivity Band 4 (9.71 um)',
    'Emissivity Band 5 (10.11 um)',
    'Emissivity Band 6 (10.58 um)',
]

LST_BAND_NAMES = ['Land Surface Temperature (K)']
QA_BAND_NAMES = ['QA Map']

# KD-tree parameters
SUBSAMPLE_FACTOR = 3
K_NEIGHBORS = 4
MAX_DISTANCE_FACTOR = 3.0
N_JOBS = -1
DEG_PER_METER = 1.0 / 111320.0


# ============================================================================
# File discovery
# ============================================================================

def discover_l2_files(input_dir):
    """
    Scan for MASTER L2 HDF5 files.
    Returns list of (filepath, granule_id) tuples.
    """
    files = []

    def scan_dir(d):
        for f in os.listdir(d):
            match = L2_PATTERN.match(f)
            if match:
                flight, line, date, start, end, ver, band, sv = match.groups()
                gid = f"{flight}_{line}_{date}_{start}_{end}_{band}_{sv}"
                files.append((os.path.join(d, f), gid))

    scan_dir(input_dir)
    for entry in os.listdir(input_dir):
        subdir = os.path.join(input_dir, entry)
        if os.path.isdir(subdir):
            scan_dir(subdir)

    return sorted(files, key=lambda x: x[1])


# ============================================================================
# HDF5 reading
# ============================================================================

def read_l2_data(filepath):
    """
    Read MASTER L2 HDF5 file.
    Returns dict with emissivity, LST, QA, and geolocation.
    """
    with h5py.File(filepath, 'r') as f:
        print(f"  Reading L2 datasets...")

        emis = f['Emissivity'][:].astype(np.float32)  # (lines, samples, bands)
        lst = f['LST'][:].astype(np.float32)           # (lines, samples)
        lat = f['Lat'][:].astype(np.float64)            # (lines, samples)
        lon = f['Lon'][:].astype(np.float64)            # (lines, samples)
        qa = f['QAmap'][:].astype(np.float32)           # (lines, samples)

    # Handle fill values
    emis[emis <= 0] = np.nan
    emis[emis > 1.5] = np.nan
    lst[lst <= 0] = np.nan
    lst[lst > 500] = np.nan

    n_lines, n_samples = lat.shape
    n_emis_bands = emis.shape[2]

    print(f"  Dimensions: {n_lines} lines x {n_samples} samples")
    print(f"  Emissivity bands: {n_emis_bands}")
    print(f"  Lat range: {lat[lat!=0].min():.4f} to {lat[lat!=0].max():.4f}")
    print(f"  Lon range: {lon[lon!=0].min():.4f} to {lon[lon!=0].max():.4f}")
    print(f"  LST range: {np.nanmin(lst):.1f} to {np.nanmax(lst):.1f} K")

    return {
        'emissivity': emis,
        'lst': lst,
        'qa': qa,
        'lat': lat,
        'lon': lon,
    }


# ============================================================================
# cKDTree reprojection
# ============================================================================

def compute_output_grid(lat, lon, resolution_m=50.0, bbox=None):
    """Compute output grid from swath lat/lon."""
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

    print(f"  Output grid: {out_width} x {out_height} pixels (~{resolution_m}m)")
    return grid_lon_2d, grid_lat_2d, transform, crs, out_width, out_height


def build_kdtree(lat, lon):
    """Build cKDTree from swath coordinates."""
    valid = np.isfinite(lat) & np.isfinite(lon) & (lat != 0) & (lon != 0)
    subsample = np.zeros_like(valid, dtype=bool)
    subsample[::SUBSAMPLE_FACTOR, ::SUBSAMPLE_FACTOR] = True
    use_mask = valid & subsample

    lat_use = lat[use_mask]
    lon_use = lon[use_mask]
    src_indices = np.where(use_mask.ravel())[0]

    print(f"  KD-tree: {len(lat_use):,} points")
    coords = np.column_stack([lat_use, lon_use])
    tree = cKDTree(coords)
    return tree, src_indices


def resample_band(tree, src_indices, grid_lat, grid_lon, swath_data, max_dist):
    """Resample a single 2D swath band onto regular grid."""
    query_points = np.column_stack([grid_lat.ravel(), grid_lon.ravel()])
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

    flat_data = swath_data.ravel()
    output = np.zeros(len(query_points), dtype=np.float32)

    for k in range(K_NEIGHBORS):
        swath_idx = src_indices[indices[:, k]]
        vals = flat_data[swath_idx]
        nan_mask = np.isnan(vals)
        vals[nan_mask] = 0
        w = weights[:, k]
        w = np.where(nan_mask, 0, w)
        output += vals * w

    output[~valid_query] = np.nan

    out_h, out_w = grid_lat.shape
    return output.reshape(out_h, out_w)


def resample_cube(tree, src_indices, grid_lat, grid_lon, swath_cube, max_dist):
    """Resample a 3D cube (lines, samples, bands) onto regular grid."""
    n_bands = swath_cube.shape[2]
    out_h, out_w = grid_lat.shape

    # Precompute query once
    query_points = np.column_stack([grid_lat.ravel(), grid_lon.ravel()])
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

    output = np.zeros((out_h * out_w, n_bands), dtype=np.float32)
    weight_accum = np.zeros((out_h * out_w, n_bands), dtype=np.float32)

    for b in range(n_bands):
        print(f"    Resampling band {b+1}/{n_bands}...")
        flat_data = swath_cube[:, :, b].ravel()
        band_out = np.zeros(len(query_points), dtype=np.float32)
        band_w = np.zeros(len(query_points), dtype=np.float32)

        for k in range(K_NEIGHBORS):
            swath_idx = src_indices[indices[:, k]]
            vals = flat_data[swath_idx]
            nan_mask = np.isnan(vals)
            vals_clean = np.where(nan_mask, 0, vals)
            w = np.where(nan_mask, 0, weights[:, k])
            band_out += vals_clean * w
            band_w += w

        band_w[band_w == 0] = np.nan
        band_out = band_out / band_w
        band_out[~valid_query] = np.nan
        output[:, b] = band_out

    return output.reshape(out_h, out_w, n_bands)


# ============================================================================
# ENVI header writing
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
            x_size = abs(transform.a)
            y_size = abs(transform.e)
            ul_x = transform.c
            ul_y = transform.f
            f.write(f'map info = {{Geographic Lat/Lon, 1, 1, {ul_x}, {ul_y}, '
                    f'{x_size}, {y_size}, WGS-84}}\n')

        if crs is not None:
            f.write(f'coordinate system string = {{{crs.to_wkt()}}}\n')

        if wavelengths is not None:
            f.write('wavelength units = Micrometers\n')
            wl_str = ', '.join(f'{w:.4f}' for w in wavelengths)
            f.write(f'wavelength = {{\n  {wl_str}}}\n')

        if band_names is not None:
            bn_str = ',\n  '.join(band_names)
            f.write(f'band names = {{\n  {bn_str}}}\n')

    return hdr_path


# ============================================================================
# Output writing
# ============================================================================

def write_output(data, output_path, crs, transform, width, height,
                 band_names, wavelengths=None, description='',
                 write_geotiff=True, write_envi=True,
                      overwrite=True):
    """Write data to GeoTIFF and/or ENVI."""
    output_files = []

    # Ensure (bands, height, width)
    if data.ndim == 2:
        stack = data[np.newaxis, :, :].astype(np.float32)
    elif data.ndim == 3 and data.shape[2] < data.shape[0]:
        stack = np.transpose(data, (2, 0, 1)).astype(np.float32)
    else:
        stack = data.astype(np.float32)

    num_bands = stack.shape[0]

    if write_geotiff:
        gtiff_path = output_path + '.tif'
        profile = {
            'driver': 'GTiff', 'dtype': 'float32',
            'width': width, 'height': height, 'count': num_bands,
            'crs': crs, 'transform': transform, 'nodata': np.nan,
            'compress': 'deflate', 'predictor': 2, 'zlevel': 6,
            'tiled': True, 'blockxsize': 256, 'blockysize': 256,
        }
        print(f"    Writing GeoTIFF: {os.path.basename(gtiff_path)}")
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
        print(f"    Writing ENVI: {os.path.basename(envi_dat)}")
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

def process_granule(filepath, granule_id, output_dir, resolution_m=50.0,
                    bbox=None,
                    
                    write_geotiff=True, write_envi=True):
    """Process a single MASTER L2 granule."""
    all_outputs = []

    data = read_l2_data(filepath)

    # Compute grid and KD-tree
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

    # --- Emissivity ---
    print(f"\n  === Processing Emissivity ===")
    emis_ortho = resample_cube(
        tree, src_indices, grid_lat, grid_lon, data['emissivity'], max_dist
    )

    n_emis = emis_ortho.shape[2]
    emis_wl = EMIS_WAVELENGTHS_UM[:n_emis]
    emis_names = EMIS_BAND_NAMES[:n_emis]

    emis_base = os.path.join(output_dir, f"MASTER_L2_EMIS_{granule_id}")
    outputs = write_output(
        emis_ortho, emis_base, crs, transform, out_w, out_h,
        band_names=emis_names, wavelengths=emis_wl,
        description=f'MASTER L2 Emissivity ({granule_id})',
        write_geotiff=write_geotiff, write_envi=write_envi,
            overwrite=overwrite
    )
    all_outputs.extend(outputs)
    del emis_ortho

    # --- LST ---
    print(f"\n  === Processing LST ===")
    lst_ortho = resample_band(
        tree, src_indices, grid_lat, grid_lon, data['lst'], max_dist
    )

    lst_base = os.path.join(output_dir, f"MASTER_L2_LST_{granule_id}")
    outputs = write_output(
        lst_ortho, lst_base, crs, transform, out_w, out_h,
        band_names=LST_BAND_NAMES,
        description=f'MASTER L2 Land Surface Temperature (K) ({granule_id})',
        write_geotiff=write_geotiff, write_envi=write_envi,
            overwrite=overwrite
    )
    all_outputs.extend(outputs)

    # --- QA ---
    print(f"\n  === Processing QA ===")
    qa_ortho = resample_band(
        tree, src_indices, grid_lat, grid_lon, data['qa'], max_dist
    )

    qa_base = os.path.join(output_dir, f"MASTER_L2_QA_{granule_id}")
    outputs = write_output(
        qa_ortho, qa_base, crs, transform, out_w, out_h,
        band_names=QA_BAND_NAMES,
        description=f'MASTER L2 QA Map ({granule_id})',
        write_geotiff=write_geotiff, write_envi=write_envi,
            overwrite=overwrite
    )
    all_outputs.extend(outputs)

    return all_outputs


def process_directory(input_dir, output_dir=None, resolution_m=50.0,
                      bbox=None,
                      
                    write_geotiff=True, write_envi=True,
                      overwrite=True):
    """Process all MASTER L2 files in a directory."""
    if output_dir is None:
        output_dir = os.path.join(input_dir, 'converted')

    print(f"Scanning for MASTER L2 files in: {input_dir}")
    files = discover_l2_files(input_dir)

    if not files:
        print("No MASTER L2 files found.")
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
                bbox=bbox,
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
        print("MASTER L2 Emissivity/LST Converter")
        print("=" * 50)
        print(f"\nUsage: python {os.path.basename(__file__)} <input_dir> [output_dir] [options]")
        print(f"\nOutput format:")
        print(f"  --geotiff         GeoTIFF only")
        print(f"  --envi            ENVI only")
        print(f"  (default: both)")
        print(f"\nResolution:")
        print(f"  --resolution <m>  Output pixel size in meters (default: 50)")
        print(f"\nExamples:")
        print(f"  python {os.path.basename(__file__)} D:\\Data\\MASTER")
        print(f"  python {os.path.basename(__file__)} D:\\Data\\MASTER --envi --resolution 25")
        print(f"\nOutputs per granule:")
        print(f"  1. Emissivity (6 TIR bands)")
        print(f"  2. LST (single band, Kelvin)")
        print(f"  3. QA Map (single band)")
        sys.exit(0)

    args = sys.argv[1:]
    write_geotiff = True
    write_envi = True
    resolution = 50.0
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

    print(f"Resolution: {resolution}m")
    print()

    process_directory(input_dir, output_dir,
                      resolution_m=resolution,
                      write_geotiff=write_geotiff, write_envi=write_envi)