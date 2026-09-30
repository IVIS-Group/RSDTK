"""
HyTES L2 Emissivity and LST Converter
Reads HyTES L2 HDF5 files, reprojects swath data using cKDTree, and outputs:
  1. Emissivity - TES-derived (164 bands, 8-11.5 µm window)
  2. LST - Land Surface Temperature (single band, Kelvin)
  3. PC Emissivity - Principal Component extended (256 bands, optional)

The L2 file does not contain geolocation — it must be read from the
matching L1 file in the same directory.

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

# Filename patterns
L2_PATTERN = re.compile(
    r'^(\d{8}t\d{6})_(.+)_L2_B(\d+)_V(\d+)\.hdf5$',
    re.IGNORECASE
)
L1_PATTERN = re.compile(
    r'^(\d{8}t\d{6})_(.+)_L1_B(\d+)_V(\d+)\.hdf5$',
    re.IGNORECASE
)

# KD-tree parameters
SUBSAMPLE_FACTOR = 3
K_NEIGHBORS = 4
MAX_DISTANCE_FACTOR = 3.0
N_JOBS = -1
DEG_PER_METER = 1.0 / 111320.0


# ============================================================================
# File discovery and matching
# ============================================================================

def discover_l2_files(input_dir):
    """
    Scan for HyTES L2 HDF5 files and match with L1 files for geolocation.
    Returns list of (l2_path, l1_path, granule_id) tuples.
    """
    l2_files = {}
    l1_files = {}

    def scan_dir(d):
        for f in os.listdir(d):
            l2_match = L2_PATTERN.match(f)
            if l2_match:
                timestamp, site, block, ver = l2_match.groups()
                key = f"{timestamp}_{site}_B{block}"
                l2_files[key] = os.path.join(d, f)
                continue
            l1_match = L1_PATTERN.match(f)
            if l1_match:
                timestamp, site, block, ver = l1_match.groups()
                key = f"{timestamp}_{site}_B{block}"
                l1_files[key] = os.path.join(d, f)

    scan_dir(input_dir)
    for entry in os.listdir(input_dir):
        subdir = os.path.join(input_dir, entry)
        if os.path.isdir(subdir):
            scan_dir(subdir)

    matched = []
    for key, l2_path in sorted(l2_files.items()):
        l1_path = l1_files.get(key)
        if l1_path is None:
            print(f"  WARNING: No matching L1 file for geolocation: {key}")
        matched.append((l2_path, l1_path, key))

    return matched


# ============================================================================
# Data reading
# ============================================================================

def read_geolocation(l1_path):
    """Read lat/lon from L1 file."""
    with h5py.File(l1_path, 'r') as f:
        lat = f['latitude'][:].astype(np.float64)
        lon = f['longitude'][:].astype(np.float64)
    print(f"  Geolocation from L1: {lat.shape}")
    return lat, lon


def read_l2_data(l2_path, include_pc=False):
    """
    Read HyTES L2 emissivity, LST, and wavelengths.

    Parameters:
        l2_path: path to L2 HDF5 file
        include_pc: if True, also read PC emissivity

    Returns dict with data arrays.
    """
    with h5py.File(l2_path, 'r') as f:
        print(f"  Reading L2 datasets...")

        emis = f['L2_Emissivity'][:].astype(np.float32)
        emis_wl = f['L2_Emissivity_Wavelengths'][:].astype(np.float64)
        lst = f['L2_LST'][:].astype(np.float32)

        attrs = dict(f.attrs)

        pc_emis = None
        pc_wl = None
        if include_pc and 'L2_Emissivity_PC' in f:
            print(f"  Reading PC Emissivity...")
            pc_emis = f['L2_Emissivity_PC'][:].astype(np.float32)
            pc_wl = f['L2_Emissivity_PC_Wavelengths'][:].astype(np.float64)

    n_lines, n_samples, n_emis = emis.shape
    print(f"  Dimensions: {n_lines} x {n_samples}")
    print(f"  Emissivity: {n_emis} bands ({emis_wl.min():.3f} - {emis_wl.max():.3f} µm)")
    print(f"  LST range: {np.nanmin(lst[lst>0]):.1f} - {np.nanmax(lst):.1f} K")

    if pc_emis is not None:
        print(f"  PC Emissivity: {pc_emis.shape[2]} bands ({pc_wl.min():.3f} - {pc_wl.max():.3f} µm)")

    # Handle fill values
    emis[emis <= 0] = np.nan
    emis[emis > 1.5] = np.nan
    lst[lst <= 0] = np.nan
    lst[lst > 1000] = np.nan

    if pc_emis is not None:
        pc_emis[pc_emis <= 0] = np.nan
        pc_emis[pc_emis > 1.5] = np.nan

    return {
        'emissivity': emis,
        'emis_wavelengths': emis_wl,
        'lst': lst,
        'pc_emissivity': pc_emis,
        'pc_wavelengths': pc_wl,
        'attributes': attrs,
    }


def estimate_resolution(lat, lon):
    """Estimate pixel resolution in meters from geolocation arrays."""
    valid = np.isfinite(lat) & np.isfinite(lon) & (lat != 0) & (lon != 0)
    mid_line = lat.shape[0] // 2
    valid_cols = np.where(valid[mid_line, :])[0]

    if len(valid_cols) > 1:
        dlat = np.abs(np.diff(lat[mid_line, valid_cols[:10]]))
        dlon = np.abs(np.diff(lon[mid_line, valid_cols[:10]]))
        center_lat = np.mean(lat[valid])
        dy = dlat.mean() / DEG_PER_METER
        dx = dlon.mean() / (DEG_PER_METER / np.cos(np.radians(center_lat)))
        res = (dx + dy) / 2.0
        print(f"  Estimated resolution: ~{res:.1f} m")
        return max(res, 1.0)
    else:
        return 5.0


# ============================================================================
# cKDTree reprojection
# ============================================================================

def compute_output_grid(lat, lon, resolution_m, bbox=None):
    """Compute output grid."""
    valid = np.isfinite(lat) & np.isfinite(lon) & (lat != 0) & (lon != 0)
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
    print(f"  Output grid: {out_w} x {out_h} pixels (~{resolution_m:.1f}m)")
    return grid_lon_2d, grid_lat_2d, transform, crs, out_w, out_h


def build_kdtree(lat, lon):
    """Build cKDTree."""
    valid = np.isfinite(lat) & np.isfinite(lon) & (lat != 0) & (lon != 0)
    subsample = np.zeros_like(valid, dtype=bool)
    subsample[::SUBSAMPLE_FACTOR, ::SUBSAMPLE_FACTOR] = True
    use_mask = valid & subsample

    coords = np.column_stack([lat[use_mask], lon[use_mask]])
    src_indices = np.where(use_mask.ravel())[0]
    tree = cKDTree(coords)
    print(f"  KD-tree: {len(coords):,} points")
    return tree, src_indices


def resample_band(tree, src_indices, grid_lat, grid_lon, swath_data, max_dist):
    """Resample single 2D band."""
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
    w_accum = np.zeros(len(query_points), dtype=np.float32)

    for k in range(K_NEIGHBORS):
        swath_idx = src_indices[indices[:, k]]
        vals = flat_data[swath_idx]
        nan_mask = np.isnan(vals)
        vals_clean = np.where(nan_mask, 0, vals)
        w = np.where(nan_mask, 0, weights[:, k])
        output += vals_clean * w
        w_accum += w

    w_accum[w_accum == 0] = np.nan
    output = output / w_accum
    output[~valid_query] = np.nan

    return output.reshape(grid_lat.shape)


def resample_cube(tree, src_indices, grid_lat, grid_lon, data_cube, max_dist):
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

    neighbor_swath_idx = np.zeros((K_NEIGHBORS, n_points), dtype=np.int64)
    for k in range(K_NEIGHBORS):
        neighbor_swath_idx[k] = src_indices[indices[:, k]]

    output = np.zeros((out_height, out_width, n_bands), dtype=np.float32)

    print(f"  Resampling {n_bands} bands (band-by-band)...")
    for b in range(n_bands):
        if (b + 1) % 50 == 0 or b == 0 or b == n_bands - 1:
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

def write_output(data, output_path, crs, transform, band_names,
                 wavelengths=None, description='',
                 write_geotiff=True, write_envi=True,
                      overwrite=True):
    """Write data to GeoTIFF and/or ENVI."""
    output_files = []

    if data.ndim == 2:
        stack = data[np.newaxis, :, :].astype(np.float32)
    elif data.ndim == 3 and data.shape[2] < data.shape[0]:
        stack = np.transpose(data, (2, 0, 1)).astype(np.float32)
    else:
        stack = data.astype(np.float32)

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
                wavelengths=wavelengths, band_names=band_names,
                description=description
        )
        output_files.append(envi_dat)
        output_files.append(envi_hdr)

    return output_files


# ============================================================================
# Main processing
# ============================================================================

def process_granule(l2_path, l1_path, granule_id, output_dir,
                    resolution_m=None, include_pc=False,
                    bbox=None, spectral_opts=None,
                    
                    write_geotiff=True, write_envi=True):
    """Process a single HyTES L2 granule."""
    all_outputs = []

    if l1_path is None:
        print(f"  ERROR: No matching L1 file for geolocation. Skipping.")
        return all_outputs

    # Read geolocation from L1
    lat, lon = read_geolocation(l1_path)

    # Read L2 data
    data = read_l2_data(l2_path, include_pc=include_pc)

    # --- Apply spectral filtering to emissivity if requested ---
    if spectral_opts is not None and spectral_opts.get('mode', 'default') != 'default':
        try:
            from spectral_tools import apply_spectral_selection
            wl_nm = data['wavelengths'] * 1000.0  # µm to nm
            good_wl = np.ones(len(wl_nm), dtype=np.int32)
            band_mask, _ = apply_spectral_selection(wl_nm, good_wl, spectral_opts)
            band_indices = np.where(band_mask)[0]
            n_orig = len(data['wavelengths'])
            data['emissivity'] = data['emissivity'][:, :, band_indices]
            data['wavelengths'] = data['wavelengths'][band_indices]
            if data.get('fwhm') is not None:
                data['fwhm'] = data['fwhm'][band_indices]
            if data.get('pc_emissivity') is not None:
                data['pc_emissivity'] = data['pc_emissivity'][:, :, band_indices]
            print(f"  Spectral filter: {len(band_indices)} of {n_orig} bands selected")
        except ImportError:
            print("  WARNING: spectral_tools.py not found, skipping spectral filter")
        except Exception as e:
            print(f"  WARNING: Spectral filtering failed: {e}")

    # Verify dimensions match
    if data['emissivity'].shape[:2] != lat.shape:
        print(f"  WARNING: L2 dimensions {data['emissivity'].shape[:2]} don't match "
              f"L1 geolocation {lat.shape}")

    # Estimate resolution
    if resolution_m is None:
        resolution_m = estimate_resolution(lat, lon)

    # Build grid and KD-tree
    print(f"\n  Computing output grid...")
    grid_lon, grid_lat, transform, crs, out_w, out_h = \
compute_output_grid(lat, lon, resolution_m, bbox=bbox)

    print(f"\n  Building KD-tree...")
    tree, src_indices = build_kdtree(lat, lon)

    center_lat = (lat[lat != 0].min() + lat[lat != 0].max()) / 2.0
    pixel_deg = resolution_m * DEG_PER_METER / np.cos(np.radians(center_lat))
    max_dist = pixel_deg * MAX_DISTANCE_FACTOR

    os.makedirs(output_dir, exist_ok=True)

    # --- Emissivity ---
    print(f"\n  === Processing Emissivity ===")
    emis_ortho = resample_cube(
        tree, src_indices, grid_lat, grid_lon, data['emissivity'], max_dist
    )

    emis_wl = data['emis_wavelengths']
    emis_names = [f'{w:.3f} um' for w in emis_wl]

    emis_base = os.path.join(output_dir, f"HyTES_L2_Emissivity_{granule_id}")
    outputs = write_output(
        emis_ortho, emis_base, crs, transform,
        band_names=emis_names, wavelengths=emis_wl,
        description=f'HyTES L2 TES Emissivity ({granule_id})',
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

    lst_base = os.path.join(output_dir, f"HyTES_L2_LST_{granule_id}")
    outputs = write_output(
        lst_ortho, lst_base, crs, transform,
        band_names=['Land Surface Temperature (K)'],
        description=f'HyTES L2 LST in Kelvin ({granule_id})',
        write_geotiff=write_geotiff, write_envi=write_envi,
            overwrite=overwrite
    )
    all_outputs.extend(outputs)

    # --- PC Emissivity (optional) ---
    if include_pc and data['pc_emissivity'] is not None:
        print(f"\n  === Processing PC Emissivity ===")
        pc_ortho = resample_cube(
            tree, src_indices, grid_lat, grid_lon, data['pc_emissivity'], max_dist
        )

        pc_wl = data['pc_wavelengths']
        pc_names = [f'{w:.3f} um' for w in pc_wl]

        pc_base = os.path.join(output_dir, f"HyTES_L2_PC_Emissivity_{granule_id}")
        outputs = write_output(
            pc_ortho, pc_base, crs, transform,
            band_names=pc_names, wavelengths=pc_wl,
            description=f'HyTES L2 PC Emissivity ({granule_id})',
            write_geotiff=write_geotiff, write_envi=write_envi,
            overwrite=overwrite
        )
        all_outputs.extend(outputs)
        del pc_ortho

    return all_outputs


def process_directory(input_dir, output_dir=None, resolution_m=None,
                      include_pc=False, bbox=None, spectral_opts=None,
                      
                    write_geotiff=True, write_envi=True,
                      overwrite=True):
    """Process all HyTES L2 files in a directory."""
    if output_dir is None:
        output_dir = os.path.join(input_dir, 'converted')

    print(f"Scanning for HyTES L2 files in: {input_dir}")
    files = discover_l2_files(input_dir)

    if not files:
        print("No HyTES L2 files found.")
        return

    print(f"Found {len(files)} file(s):\n")

    all_outputs = []
    for i, (l2_path, l1_path, gid) in enumerate(files, 1):
        print(f"{'='*70}")
        print(f"Granule {i}/{len(files)}: {gid}")
        print(f"  L2: {os.path.basename(l2_path)}")
        print(f"  L1: {os.path.basename(l1_path) if l1_path else 'NOT FOUND'}")
        print(f"{'='*70}")

        try:
            outputs = process_granule(
                l2_path, l1_path, gid, output_dir,
                resolution_m=resolution_m,
                include_pc=include_pc,
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
        print("HyTES L2 Emissivity/LST Converter")
        print("=" * 50)
        print(f"\nUsage: python {os.path.basename(__file__)} <input_dir> [output_dir] [options]")
        print(f"\nOutput format:")
        print(f"  --geotiff         GeoTIFF only")
        print(f"  --envi            ENVI only")
        print(f"  (default: both)")
        print(f"\nOptions:")
        print(f"  --pc-emis         Include PC emissivity (256 bands, full spectral range)")
        print(f"  --resolution <m>  Output pixel size in meters (default: auto-estimate)")
        print(f"\nThe L1 file must be in the same directory for geolocation.")
        print(f"\nOutputs per granule:")
        print(f"  1. Emissivity (TES-derived, ~164 bands, 8-11.5 µm)")
        print(f"  2. LST (single band, Kelvin)")
        print(f"  3. PC Emissivity (256 bands, optional)")
        sys.exit(0)

    args = sys.argv[1:]
    write_geotiff = True
    write_envi = True
    resolution = None
    include_pc = False
    positional = []

    i = 0
    while i < len(args):
        a = args[i].lower()
        if a == '--geotiff':
            write_geotiff = True; write_envi = False
        elif a == '--envi':
            write_geotiff = False; write_envi = True
        elif a == '--pc-emis':
            include_pc = True
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

    print(f"PC Emissivity: {'yes' if include_pc else 'no'}")
    print()

    process_directory(input_dir, output_dir, resolution_m=resolution,
                      include_pc=include_pc,
                      write_geotiff=write_geotiff, write_envi=write_envi)