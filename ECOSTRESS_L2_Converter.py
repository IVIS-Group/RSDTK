"""
ECOSTRESS L2 LSTE Swath Converter
Reads ECO_L2_LSTE (swath) HDF5 + matching L1B_GEO file, applies scale factors,
reprojects to regular grid using cKDTree, and outputs three files per granule:
  1. Emissivity: Emis1-5 (5 bands)
  2. LST: Land Surface Temperature + Wideband Emissivity (2 bands)
  3. QA/Errors: Emis1_err-5_err + LST_err + cloud_mask + water_mask (8 bands)

Supports Geographic (WGS84) and UTM output projections at 70m native resolution.
Supports subsetting via vector files, text files, manual bbox, or auto-tiling.

ECOSTRESS TIR Band Wavelengths:
  Band 1: 8.29 µm
  Band 2: 8.63 µm
  Band 3: 9.07 µm
  Band 4: 10.49 µm
  Band 5: 12.01 µm
"""

import os
import re
import sys
import glob
import numpy as np
import h5py
from scipy.spatial import cKDTree
from collections import defaultdict
import warnings

try:
    import rasterio
    from rasterio.transform import from_bounds
    from rasterio.crs import CRS
except ImportError:
    print("ERROR: rasterio is required. Install with: pip install rasterio")
    sys.exit(1)

# Import subsetting module
try:
    from subset_tool import (
        BoundingBox, bbox_from_file, bbox_from_manual,
        parse_subset_args, TileGrid, compute_tile_grid, get_swath_bbox
    )
    SUBSET_AVAILABLE = True
except ImportError:
    print("WARNING: subset_tool.py not found. Subsetting/tiling disabled.")
    print("Place subset_tool.py in the same directory as this script.")
    SUBSET_AVAILABLE = False

# ============================================================================
# Constants and configuration
# ============================================================================

# Filename patterns
LSTE_PATTERN = re.compile(
    r'^ECOv002_L2_LSTE_(\d+)_(\d+)_(\d{8}T\d{6})_(\d{4})_(\d{2})\.h5$',
    re.IGNORECASE
)
GEO_PATTERN = re.compile(
    r'^ECOv002_L1B_GEO_(\d+)_(\d+)_(\d{8}T\d{6})_(\d{4})_(\d{2})\.h5$',
    re.IGNORECASE
)

# HDF5 dataset paths
GEO_LAT_PATH = 'Geolocation/latitude'
GEO_LON_PATH = 'Geolocation/longitude'

SDS_PREFIX = 'SDS/'

# Emissivity datasets and scale factors (from User Guide Table 3)
EMIS_BANDS = {
    'Emis1': {'path': f'{SDS_PREFIX}Emis1', 'scale': 0.002, 'offset': 0.49, 'fill': 0},
    'Emis2': {'path': f'{SDS_PREFIX}Emis2', 'scale': 0.002, 'offset': 0.49, 'fill': 0},
    'Emis3': {'path': f'{SDS_PREFIX}Emis3', 'scale': 0.002, 'offset': 0.49, 'fill': 0},
    'Emis4': {'path': f'{SDS_PREFIX}Emis4', 'scale': 0.002, 'offset': 0.49, 'fill': 0},
    'Emis5': {'path': f'{SDS_PREFIX}Emis5', 'scale': 0.002, 'offset': 0.49, 'fill': 0},
}

LST_BANDS = {
    'LST': {'path': f'{SDS_PREFIX}LST', 'scale': 0.02, 'offset': 0.0, 'fill': 0},
    'EmisWB': {'path': f'{SDS_PREFIX}EmisWB', 'scale': 0.002, 'offset': 0.49, 'fill': 0},
}

QA_BANDS = {
    'Emis1_err': {'path': f'{SDS_PREFIX}Emis1_err', 'scale': 0.0001, 'offset': 0.0, 'fill': 0},
    'Emis2_err': {'path': f'{SDS_PREFIX}Emis2_err', 'scale': 0.0001, 'offset': 0.0, 'fill': 0},
    'Emis3_err': {'path': f'{SDS_PREFIX}Emis3_err', 'scale': 0.0001, 'offset': 0.0, 'fill': 0},
    'Emis4_err': {'path': f'{SDS_PREFIX}Emis4_err', 'scale': 0.0001, 'offset': 0.0, 'fill': 0},
    'Emis5_err': {'path': f'{SDS_PREFIX}Emis5_err', 'scale': 0.0001, 'offset': 0.0, 'fill': 0},
    'LST_err':   {'path': f'{SDS_PREFIX}LST_err',   'scale': 0.04,   'offset': 0.0, 'fill': 0},
    'cloud_mask': {'path': f'{SDS_PREFIX}cloud_mask', 'scale': 1, 'offset': 0, 'fill': 255},
    'water_mask': {'path': f'{SDS_PREFIX}water_mask', 'scale': 1, 'offset': 0, 'fill': 255},
}

# ECOSTRESS TIR wavelengths in micrometers
EMIS_WAVELENGTHS_UM = [8.29, 8.63, 9.07, 10.49, 12.01]

# Band names for ENVI headers
EMIS_BAND_NAMES = [
    'Emissivity Band 1 (8.29 um)',
    'Emissivity Band 2 (8.63 um)',
    'Emissivity Band 3 (9.07 um)',
    'Emissivity Band 4 (10.49 um)',
    'Emissivity Band 5 (12.01 um)',
]

LST_BAND_NAMES = ['Land Surface Temperature (K)', 'Wideband Emissivity']

QA_BAND_NAMES = [
    'Emissivity Error Band 1 (8.29 um)',
    'Emissivity Error Band 2 (8.63 um)',
    'Emissivity Error Band 3 (9.07 um)',
    'Emissivity Error Band 4 (10.49 um)',
    'Emissivity Error Band 5 (12.01 um)',
    'LST Error (K)',
    'Cloud Mask',
    'Water Mask',
]

# Native resolution in meters
NATIVE_RESOLUTION_M = 70.0

# Approximate degrees per meter at equator
DEG_PER_METER_EQUATOR = 1.0 / 111320.0

# cKDTree resampling parameters
SUBSAMPLE_FACTOR = 3
EDGE_PRESERVE_COLS = 270
K_NEIGHBORS = 4
MAX_DISTANCE_FACTOR = 3.0
N_JOBS = -1


# ============================================================================
# File discovery and matching
# ============================================================================

def discover_files(input_dir):
    """
    Scan input directory for L2_LSTE and L1B_GEO HDF5 files.
    Match them by orbit_scene_timestamp.
    Returns list of (lste_path, geo_path, granule_id) tuples.
    """
    lste_files = {}
    geo_files = {}

    for filename in os.listdir(input_dir):
        if not filename.lower().endswith('.h5'):
            continue

        lste_match = LSTE_PATTERN.match(filename)
        if lste_match:
            orbit, scene, timestamp, build, iteration = lste_match.groups()
            key = f"{orbit}_{scene}_{timestamp}"
            lste_files[key] = os.path.join(input_dir, filename)
            continue

        geo_match = GEO_PATTERN.match(filename)
        if geo_match:
            orbit, scene, timestamp, build, iteration = geo_match.groups()
            key = f"{orbit}_{scene}_{timestamp}"
            geo_files[key] = os.path.join(input_dir, filename)

    matched = []
    for key, lste_path in sorted(lste_files.items()):
        if key in geo_files:
            matched.append((lste_path, geo_files[key], key))
        else:
            print(f"  WARNING: No matching GEO file for LSTE granule: {key}")

    return matched


# ============================================================================
# HDF5 reading
# ============================================================================

def read_geo_data(geo_path):
    """Read latitude and longitude arrays from GEO file."""
    print(f"  Reading geolocation: {os.path.basename(geo_path)}")
    with h5py.File(geo_path, 'r') as f:
        lat = f[GEO_LAT_PATH][:]
        lon = f[GEO_LON_PATH][:]
    print(f"  Geolocation shape: {lat.shape}")
    return lat, lon


def read_sds_band(h5file, band_info):
    """
    Read a single SDS band from an open HDF5 file and apply scale/offset.
    Returns float32 array with fill values set to NaN.
    """
    ds = h5file[band_info['path']]
    raw = ds[:]

    data = raw.astype(np.float32)
    fill_val = band_info['fill']
    fill_mask = (raw == fill_val)

    if band_info['scale'] != 1 or band_info['offset'] != 0:
        data = data * band_info['scale'] + band_info['offset']

    data[fill_mask] = np.nan
    return data


def read_band_group(lste_path, band_dict):
    """
    Read a group of bands from the LSTE file.
    Returns list of (band_name, data_array) tuples.
    """
    bands = []
    with h5py.File(lste_path, 'r') as f:
        for band_name, band_info in band_dict.items():
            print(f"    Reading {band_name}...")
            data = read_sds_band(f, band_info)
            bands.append((band_name, data))
    return bands


# ============================================================================
# cKDTree reprojection
# ============================================================================

def compute_output_grid(lat, lon, projection='geographic', subset_bbox=None):
    """
    Compute the output grid parameters based on the swath extent,
    optionally clipped to a subset bounding box.

    Parameters:
        lat, lon: 2D arrays of swath coordinates
        projection: 'geographic' or 'utm'
        subset_bbox: BoundingBox to clip to (optional)

    Returns:
        grid_lon, grid_lat: 2D arrays of output grid coordinates
        out_crs: rasterio CRS object
        transform: rasterio Affine transform
        out_width, out_height: output dimensions
    """
    valid = np.isfinite(lat) & np.isfinite(lon) & (lat != 0) & (lon != 0)

    if not valid.any():
        raise ValueError("No valid geolocation data found.")

    lat_valid = lat[valid]
    lon_valid = lon[valid]

    lat_min, lat_max = lat_valid.min(), lat_valid.max()
    lon_min, lon_max = lon_valid.min(), lon_valid.max()

    # Apply subset bounding box if provided
    if subset_bbox is not None:
        swath_bbox = BoundingBox(lat_min, lat_max, lon_min, lon_max)
        clipped = swath_bbox.intersection(subset_bbox)
        if clipped is None:
            raise ValueError("Subset bounding box does not overlap with swath extent.\n"
                             f"  Swath:  {swath_bbox}\n"
                             f"  Subset: {subset_bbox}")
        lat_min, lat_max = clipped.min_lat, clipped.max_lat
        lon_min, lon_max = clipped.min_lon, clipped.max_lon
        print(f"  Clipped extent: Lat [{lat_min:.4f}, {lat_max:.4f}], "
              f"Lon [{lon_min:.4f}, {lon_max:.4f}]")
    else:
        print(f"  Swath extent: Lat [{lat_min:.4f}, {lat_max:.4f}], "
              f"Lon [{lon_min:.4f}, {lon_max:.4f}]")

    if projection == 'utm':
        center_lon = (lon_min + lon_max) / 2.0
        center_lat = (lat_min + lat_max) / 2.0
        utm_zone = int((center_lon + 180) / 6) + 1
        is_south = center_lat < 0
        epsg = 32600 + utm_zone if not is_south else 32700 + utm_zone
        out_crs = CRS.from_epsg(epsg)

        print(f"  UTM Zone: {utm_zone}{'S' if is_south else 'N'} (EPSG:{epsg})")

        from pyproj import Transformer
        transformer = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)

        # Transform corners to get UTM extent
        corners_lon = [lon_min, lon_min, lon_max, lon_max]
        corners_lat = [lat_min, lat_max, lat_min, lat_max]
        x_corners, y_corners = transformer.transform(corners_lon, corners_lat)
        x_min, x_max = min(x_corners), max(x_corners)
        y_min, y_max = min(y_corners), max(y_corners)

        res_x = NATIVE_RESOLUTION_M
        res_y = NATIVE_RESOLUTION_M

        grid_x = np.arange(x_min, x_max, res_x)
        grid_y = np.arange(y_max, y_min, -res_y)

        out_width = len(grid_x)
        out_height = len(grid_y)

        transform = rasterio.transform.from_origin(x_min, y_max, res_x, res_y)

        grid_xx, grid_yy = np.meshgrid(grid_x, grid_y)

        transformer_inv = Transformer.from_crs(f"EPSG:{epsg}", "EPSG:4326", always_xy=True)
        grid_lon_2d, grid_lat_2d = transformer_inv.transform(grid_xx, grid_yy)

    else:  # geographic
        out_crs = CRS.from_epsg(4326)

        center_lat = (lat_min + lat_max) / 2.0
        res_y_deg = NATIVE_RESOLUTION_M * DEG_PER_METER_EQUATOR
        res_x_deg = NATIVE_RESOLUTION_M * DEG_PER_METER_EQUATOR / np.cos(np.radians(center_lat))

        print(f"  Output resolution: {res_x_deg:.6f} x {res_y_deg:.6f} degrees "
              f"(~{NATIVE_RESOLUTION_M}m)")

        grid_lon_1d = np.arange(lon_min, lon_max, res_x_deg)
        grid_lat_1d = np.arange(lat_max, lat_min, -res_y_deg)

        out_width = len(grid_lon_1d)
        out_height = len(grid_lat_1d)

        transform = rasterio.transform.from_origin(lon_min, lat_max, res_x_deg, res_y_deg)

        grid_lon_2d, grid_lat_2d = np.meshgrid(grid_lon_1d, grid_lat_1d)

    print(f"  Output grid: {out_width} x {out_height} pixels")

    return grid_lon_2d, grid_lat_2d, out_crs, transform, out_width, out_height


def build_kdtree(lat, lon):
    """
    Build a cKDTree from swath lat/lon with subsampling for performance.
    Preserves edge pixels to avoid gaps from ECOSTRESS scan geometry.
    """
    nrows, ncols = lat.shape
    print(f"  Building KD-tree from swath ({nrows} x {ncols})...")

    valid = np.isfinite(lat) & np.isfinite(lon) & (lat != 0) & (lon != 0)

    subsample_mask = np.zeros_like(valid, dtype=bool)
    subsample_mask[::SUBSAMPLE_FACTOR, ::SUBSAMPLE_FACTOR] = True

    edge_mask = np.zeros_like(valid, dtype=bool)
    edge_mask[:, :EDGE_PRESERVE_COLS] = True
    edge_mask[:, -EDGE_PRESERVE_COLS:] = True

    use_mask = valid & (subsample_mask | edge_mask)

    lat_use = lat[use_mask]
    lon_use = lon[use_mask]
    src_indices = np.where(use_mask.ravel())[0]

    print(f"  KD-tree points: {len(lat_use):,} "
          f"(from {valid.sum():,} valid, subsample={SUBSAMPLE_FACTOR}x, "
          f"edge={EDGE_PRESERVE_COLS} cols preserved)")

    coords = np.column_stack([lat_use, lon_use])
    tree = cKDTree(coords)

    return tree, src_indices, use_mask


def resample_bands(tree, src_indices, lat, lon, grid_lat, grid_lon,
                   band_data_list, max_dist_deg):
    """
    Resample multiple bands from swath to regular grid using prebuilt KD-tree.
    """
    out_height, out_width = grid_lat.shape

    query_points = np.column_stack([grid_lat.ravel(), grid_lon.ravel()])

    print(f"  Querying KD-tree ({len(query_points):,} grid points, k={K_NEIGHBORS})...")
    distances, indices = tree.query(
        query_points,
        k=K_NEIGHBORS,
        workers=N_JOBS
    )

    if K_NEIGHBORS == 1:
        distances = distances[:, np.newaxis]
        indices = indices[:, np.newaxis]

    valid_query = distances[:, 0] < max_dist_deg

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        weights = 1.0 / distances
        weights[~np.isfinite(weights)] = 1e10

    weight_sum = weights.sum(axis=1, keepdims=True)
    weight_sum[weight_sum == 0] = 1.0
    weights = weights / weight_sum

    resampled = []
    for band_name, swath_data in band_data_list:
        print(f"    Resampling {band_name}...")
        flat_swath = swath_data.ravel()

        neighbor_vals = np.full_like(distances, np.nan, dtype=np.float32)
        for k in range(K_NEIGHBORS):
            swath_flat_idx = src_indices[indices[:, k]]
            neighbor_vals[:, k] = flat_swath[swath_flat_idx]

        nan_mask = np.isnan(neighbor_vals)
        masked_weights = weights.copy()
        masked_weights[nan_mask] = 0.0

        w_sum = masked_weights.sum(axis=1)
        w_sum[w_sum == 0] = np.nan

        weighted_vals = np.nansum(neighbor_vals * masked_weights, axis=1) / w_sum
        weighted_vals[~valid_query] = np.nan

        resampled.append((band_name, weighted_vals.reshape(out_height, out_width)))

    return resampled


# ============================================================================
# ENVI header writing
# ============================================================================

def build_envi_header(envi_path, height, width, num_bands, dtype, crs, transform,
                      wavelengths=None, band_names=None, description='',
                      interleave='bsq'):
    """Write a custom ENVI header file with wavelength metadata."""
    dtype_map = {
        'uint8': 1, 'int16': 2, 'int32': 3, 'float32': 4,
        'float64': 5, 'uint16': 12, 'uint32': 13,
    }

    dtype_str = np.dtype(dtype).name
    envi_dtype = dtype_map.get(dtype_str, 4)

    hdr_path = envi_path + '.hdr' if not envi_path.endswith('.hdr') else envi_path

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

            if crs and crs.is_geographic:
                f.write(f'map info = {{Geographic Lat/Lon, 1, 1, {ul_x}, {ul_y}, '
                        f'{x_size}, {y_size}, WGS-84}}\n')
            elif crs and crs.is_projected:
                f.write(f'map info = {{UTM, 1, 1, {ul_x}, {ul_y}, '
                        f'{x_size}, {y_size}, WGS-84}}\n')

        if crs is not None:
            f.write(f'coordinate system string = {{{crs.to_wkt()}}}\n')

        if wavelengths is not None and len(wavelengths) > 0:
            f.write('wavelength units = Micrometers\n')
            wl_all = list(wavelengths)
            while len(wl_all) < num_bands:
                wl_all.append(0.00)
            wl_str = ', '.join(f'{w:.2f}' for w in wl_all)
            f.write(f'wavelength = {{\n  {wl_str}}}\n')

        if band_names is not None:
            bn_str = ',\n  '.join(band_names)
            f.write(f'band names = {{\n  {bn_str}}}\n')

    return hdr_path


# ============================================================================
# Output writing
# ============================================================================

def write_output(data_bands, output_path, crs, transform, width, height,
                 band_names, wavelengths=None, description='',
                 write_geotiff=True, write_envi=True):
    """
    Write resampled bands to GeoTIFF and/or ENVI format.
    """
    output_files = []
    num_bands = len(data_bands)

    stack = np.zeros((num_bands, height, width), dtype=np.float32)
    for i, (bname, bdata) in enumerate(data_bands):
        stack[i] = bdata

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

        print(f"    Writing GeoTIFF: {os.path.basename(gtiff_path)}")
        with rasterio.open(gtiff_path, 'w', **profile) as dst:
            for i in range(num_bands):
                dst.write(stack[i], i + 1)
                dst.set_band_description(i + 1, band_names[i])
        output_files.append(gtiff_path)

    if write_envi:
        envi_dat_path = output_path + '.dat'
        envi_hdr_path = output_path + '.hdr'

        print(f"    Writing ENVI: {os.path.basename(envi_dat_path)}")
        stack.tofile(envi_dat_path)

        build_envi_header(
            envi_hdr_path,
            height=height,
            width=width,
            num_bands=num_bands,
            dtype=np.float32,
            crs=crs,
            transform=transform,
            wavelengths=wavelengths,
            band_names=band_names,
            description=description,
            interleave='bsq'
        )
        output_files.append(envi_dat_path)
        output_files.append(envi_hdr_path)

    return output_files


# ============================================================================
# Main processing
# ============================================================================

def process_single_extent(lste_path, geo_path, granule_id, output_dir,
                          lat, lon, tree, src_indices,
                          projection='geographic', subset_bbox=None,
                          tile_label=None,
                          write_geotiff=True, write_envi=True):
    """
    Process a single extent (full swath, subset, or tile) for one granule.
    The KD-tree is prebuilt and shared across tiles/subsets.

    Parameters:
        lste_path, geo_path: paths to data files
        granule_id: string identifier
        output_dir: output directory
        lat, lon: full swath coordinate arrays
        tree, src_indices: prebuilt KD-tree and index mapping
        projection: 'geographic' or 'utm'
        subset_bbox: BoundingBox to clip to (optional)
        tile_label: string label for tile (e.g., 'R0C1') or None
        write_geotiff, write_envi: output format flags

    Returns:
        list of output file paths
    """
    all_outputs = []

    # File name suffix for tiles
    suffix = f"_{tile_label}" if tile_label else ""

    # Compute output grid (clipped if subset_bbox provided)
    grid_lon, grid_lat, out_crs, transform, out_width, out_height = \
        compute_output_grid(lat, lon, projection=projection, subset_bbox=subset_bbox)

    if out_width == 0 or out_height == 0:
        print(f"  WARNING: Empty grid for extent{suffix}, skipping.")
        return all_outputs

    # Compute max search distance
    if projection == 'geographic':
        center_lat = (lat[np.isfinite(lat)].min() + lat[np.isfinite(lat)].max()) / 2.0
        pixel_deg = NATIVE_RESOLUTION_M * DEG_PER_METER_EQUATOR / np.cos(np.radians(center_lat))
    else:
        pixel_deg = NATIVE_RESOLUTION_M * DEG_PER_METER_EQUATOR
    max_dist = pixel_deg * MAX_DISTANCE_FACTOR

    # --- Emissivity ---
    print(f"\n  === Processing Emissivity{suffix} ===")
    emis_data = read_band_group(lste_path, EMIS_BANDS)
    emis_resampled = resample_bands(
        tree, src_indices, lat, lon, grid_lat, grid_lon, emis_data, max_dist
    )

    emis_base = os.path.join(output_dir, f"ECO_L2_EMIS_{granule_id}{suffix}")
    outputs = write_output(
        emis_resampled, emis_base, out_crs, transform, out_width, out_height,
        band_names=EMIS_BAND_NAMES,
        wavelengths=EMIS_WAVELENGTHS_UM,
        description='ECOSTRESS L2 Emissivity (Bands 1-5)',
        write_geotiff=write_geotiff,
        write_envi=write_envi
    )
    all_outputs.extend(outputs)

    # --- LST ---
    print(f"\n  === Processing LST{suffix} ===")
    lst_data = read_band_group(lste_path, LST_BANDS)
    lst_resampled = resample_bands(
        tree, src_indices, lat, lon, grid_lat, grid_lon, lst_data, max_dist
    )

    lst_base = os.path.join(output_dir, f"ECO_L2_LST_{granule_id}{suffix}")
    outputs = write_output(
        lst_resampled, lst_base, out_crs, transform, out_width, out_height,
        band_names=LST_BAND_NAMES,
        wavelengths=None,
        description='ECOSTRESS L2 Land Surface Temperature (K) + Wideband Emissivity',
        write_geotiff=write_geotiff,
        write_envi=write_envi
    )
    all_outputs.extend(outputs)

    # --- QA/Errors ---
    print(f"\n  === Processing QA/Errors{suffix} ===")
    qa_data = read_band_group(lste_path, QA_BANDS)
    qa_resampled = resample_bands(
        tree, src_indices, lat, lon, grid_lat, grid_lon, qa_data, max_dist
    )

    qa_base = os.path.join(output_dir, f"ECO_L2_QA_{granule_id}{suffix}")
    outputs = write_output(
        qa_resampled, qa_base, out_crs, transform, out_width, out_height,
        band_names=QA_BAND_NAMES,
        wavelengths=None,
        description='ECOSTRESS L2 QA/Error Bands',
        write_geotiff=write_geotiff,
        write_envi=write_envi
    )
    all_outputs.extend(outputs)

    return all_outputs


def process_granule(lste_path, geo_path, granule_id, output_dir,
                    projection='geographic', subset_opts=None,
                    write_geotiff=True, write_envi=True):
    """
    Process a single ECOSTRESS L2 LSTE granule with optional subsetting/tiling.

    Parameters:
        lste_path: path to L2_LSTE HDF5 file
        geo_path: path to L1B_GEO HDF5 file
        granule_id: string identifier
        output_dir: output directory
        projection: 'geographic' or 'utm'
        subset_opts: dict from parse_subset_args() or None for full extent
        write_geotiff, write_envi: output format flags

    Returns:
        list of output file paths
    """
    all_outputs = []

    # Read geolocation (always need full swath for KD-tree)
    lat, lon = read_geo_data(geo_path)

    # Build KD-tree once (reused across all tiles/subsets)
    tree, src_indices, use_mask = build_kdtree(lat, lon)

    # Determine processing mode
    if subset_opts is None or subset_opts['mode'] == 'full':
        # Full extent
        print(f"\n  Mode: Full extent")
        outputs = process_single_extent(
            lste_path, geo_path, granule_id, output_dir,
            lat, lon, tree, src_indices,
            projection=projection,
            write_geotiff=write_geotiff,
            write_envi=write_envi
        )
        all_outputs.extend(outputs)

    elif subset_opts['mode'] == 'subset':
        # Bounding box subset
        bbox = subset_opts['bbox']
        print(f"\n  Mode: Subset to {bbox}")
        outputs = process_single_extent(
            lste_path, geo_path, granule_id, output_dir,
            lat, lon, tree, src_indices,
            projection=projection,
            subset_bbox=bbox,
            tile_label='subset',
            write_geotiff=write_geotiff,
            write_envi=write_envi
        )
        all_outputs.extend(outputs)

    elif subset_opts['mode'] == 'tile':
        # Auto-tiling
        tile_size = subset_opts.get('tile_size', 2048)
        print(f"\n  Mode: Auto-tiling ({tile_size}x{tile_size} pixels)")

        swath_bbox = get_swath_bbox(lat, lon)
        tile_grid = TileGrid(swath_bbox, tile_size_pixels=tile_size,
                             resolution_m=NATIVE_RESOLUTION_M)
        print(f"  {tile_grid}")

        for row, col, tile_bbox in tile_grid.get_all_tiles():
            tile_label = f"R{row}C{col}"
            print(f"\n  --- Tile {tile_label} ---")
            print(f"  {tile_bbox}")

            try:
                outputs = process_single_extent(
                    lste_path, geo_path, granule_id, output_dir,
                    lat, lon, tree, src_indices,
                    projection=projection,
                    subset_bbox=tile_bbox,
                    tile_label=tile_label,
                    write_geotiff=write_geotiff,
                    write_envi=write_envi
                )
                all_outputs.extend(outputs)
            except ValueError as e:
                print(f"  Skipping tile {tile_label}: {e}")

    return all_outputs


def process_directory(input_dir, output_dir=None, projection='geographic',
                      subset_opts=None, write_geotiff=True, write_envi=True,
                      overwrite=True):
    """
    Process all matched L2_LSTE + L1B_GEO granule pairs in a directory.
    """
    if output_dir is None:
        output_dir = os.path.join(input_dir, 'converted')

    os.makedirs(output_dir, exist_ok=True)

    print(f"Scanning for ECOSTRESS L2 LSTE files in: {input_dir}")
    matched = discover_files(input_dir)

    if not matched:
        print("No matched LSTE + GEO file pairs found.")
        print("Ensure both ECOv002_L2_LSTE_*.h5 and ECOv002_L1B_GEO_*.h5 files "
              "are in the same directory with matching orbit/scene/timestamp.")
        return

    print(f"Found {len(matched)} matched granule pair(s):\n")

    all_outputs = []
    for i, (lste_path, geo_path, granule_id) in enumerate(matched, 1):
        print(f"{'='*70}")
        print(f"Granule {i}/{len(matched)}: {granule_id}")
        print(f"  LSTE: {os.path.basename(lste_path)}")
        print(f"  GEO:  {os.path.basename(geo_path)}")
        print(f"{'='*70}")

        try:
            outputs = process_granule(
                lste_path, geo_path, granule_id, output_dir,
                projection=projection,
                subset_opts=subset_opts,
                write_geotiff=write_geotiff,
                write_envi=write_envi
            )
            all_outputs.extend(outputs)
        except Exception as e:
            print(f"  ERROR processing granule {granule_id}: {e}")
            import traceback
            traceback.print_exc()

        print()

    print(f"Processing complete. {len(all_outputs)} output file(s) created in: {output_dir}")


# ============================================================================
# Command-line interface
# ============================================================================

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("ECOSTRESS L2 LSTE Swath Converter")
        print("=" * 50)
        print(f"\nUsage: python {os.path.basename(__file__)} <input_dir> [output_dir] [options]")
        print(f"\nProjection options:")
        print(f"  --geographic       Output in Geographic WGS84 (default)")
        print(f"  --utm              Output in UTM projection")
        print(f"\nOutput format options:")
        print(f"  --geotiff          Output GeoTIFF only")
        print(f"  --envi             Output ENVI only")
        print(f"  (default: both GeoTIFF and ENVI)")
        print(f"\nSubsetting options:")
        print(f"  --subset <file>    Clip to extent of vector/text file")
        print(f"                     Supported: .shp, .kml, .kmz, .geojson, .json, .txt, .csv")
        print(f"  --bbox <min_lat> <max_lat> <min_lon> <max_lon>")
        print(f"                     Clip to manual bounding box")
        print(f"  --tile [size]      Auto-tile output (default: 2048x2048 pixels)")
        print(f"  (default: full swath extent)")
        print(f"\nText file bbox format:")
        print(f"  min_lat = -15.5")
        print(f"  max_lat = -14.0")
        print(f"  min_lon = 165.0")
        print(f"  max_lon = 167.5")
        print(f"\nExamples:")
        print(f"  python {os.path.basename(__file__)} D:\\Data\\ECOSTRESS")
        print(f"  python {os.path.basename(__file__)} D:\\Data\\ECOSTRESS --utm --geotiff")
        print(f"  python {os.path.basename(__file__)} D:\\Data\\ECOSTRESS --subset study_area.shp")
        print(f"  python {os.path.basename(__file__)} D:\\Data\\ECOSTRESS --bbox -15.5 -14.0 165.0 167.5")
        print(f"  python {os.path.basename(__file__)} D:\\Data\\ECOSTRESS --tile 2048")
        print(f"\nOutputs per granule (per tile/subset):")
        print(f"  1. Emissivity (5 bands): Emis1-5")
        print(f"  2. LST (2 bands): LST + Wideband Emissivity")
        print(f"  3. QA/Errors (8 bands): Emis1-5_err + LST_err + cloud + water masks")
        sys.exit(0)

    # Parse arguments
    args = sys.argv[1:]
    projection = 'geographic'
    write_geotiff = True
    write_envi = True
    positional = []
    subset_opts = None

    # Extract projection and format flags first
    filtered_args = []
    i = 0
    while i < len(args):
        arg = args[i].lower()
        if arg == '--utm':
            projection = 'utm'
            i += 1
        elif arg == '--geographic':
            projection = 'geographic'
            i += 1
        elif arg == '--geotiff':
            write_geotiff = True
            write_envi = False
            i += 1
        elif arg == '--envi':
            write_geotiff = False
            write_envi = True
            i += 1
        else:
            filtered_args.append(args[i])
            i += 1

    # Parse subsetting options from remaining args
    if SUBSET_AVAILABLE:
        subset_opts = parse_subset_args(filtered_args)
        # Remove subset-related args from positional
        skip_next = 0
        for j, arg in enumerate(filtered_args):
            if skip_next > 0:
                skip_next -= 1
                continue
            a = arg.lower()
            if a == '--subset':
                skip_next = 1
            elif a == '--bbox':
                skip_next = 4
            elif a == '--tile':
                if j + 1 < len(filtered_args) and filtered_args[j + 1].isdigit():
                    skip_next = 1
            elif not a.startswith('--'):
                positional.append(arg)
    else:
        positional = [a for a in filtered_args if not a.startswith('--')]

    if not positional:
        print("Error: No input directory specified.")
        sys.exit(1)

    input_dir = positional[0]
    output_dir = positional[1] if len(positional) > 1 else None

    if not os.path.isdir(input_dir):
        print(f"Error: Input directory does not exist: {input_dir}")
        sys.exit(1)

    print(f"Projection: {projection}")
    print(f"Output formats: {'GeoTIFF' if write_geotiff else ''}"
          f"{' + ' if write_geotiff and write_envi else ''}"
          f"{'ENVI' if write_envi else ''}")
    if subset_opts:
        print(f"Subsetting mode: {subset_opts['mode']}")
    print()

    process_directory(
        input_dir, output_dir,
        projection=projection,
        subset_opts=subset_opts,
        write_geotiff=write_geotiff,
        write_envi=write_envi
    )