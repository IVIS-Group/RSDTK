"""
VIIRS Data Converter v1.0
Converts NASA VIIRS swath products (VNP02IMG, VNP02MOD, VNP21) to GeoTIFF and ENVI formats.
Handles swath-to-grid reprojection with user-selectable projection and resampling.

Supported Products:
  - VNP02IMG / VJ102IMG: L1B Imagery Resolution Bands (I01-I05), 375m
  - VNP02MOD / VJ102MOD: L1B Moderate Resolution Bands (M01-M16), 750m
  - VNP21 / VJ121: Land Surface Temperature & Emissivity, 750m

Requirements:
  pip install customtkinter numpy netCDF4 rasterio pyproj scipy
"""

import os
import re
import sys
import glob
import threading
import traceback
import numpy as np
from datetime import datetime

import customtkinter as ctk
from tkinter import filedialog, messagebox

try:
    import netCDF4 as nc4
except ImportError:
    nc4 = None

try:
    import rasterio
    from rasterio.transform import from_bounds
    from rasterio.crs import CRS
except ImportError:
    rasterio = None

try:
    from pyproj import Transformer, CRS as PyprojCRS
except ImportError:
    Transformer = None

try:
    from scipy.interpolate import griddata
    from scipy.spatial import cKDTree
except ImportError:
    griddata = None
    cKDTree = None

# =============================================================================
# VIIRS Product Definitions
# =============================================================================

VIIRS_BAND_INFO = {
    "I01": {"wavelength": 0.640, "type": "RSB", "resolution": 375},
    "I02": {"wavelength": 0.865, "type": "RSB", "resolution": 375},
    "I03": {"wavelength": 1.610, "type": "RSB", "resolution": 375},
    "I04": {"wavelength": 3.740, "type": "TEB", "resolution": 375},
    "I05": {"wavelength": 11.450, "type": "TEB", "resolution": 375},
    "M01": {"wavelength": 0.412, "type": "RSB", "resolution": 750},
    "M02": {"wavelength": 0.445, "type": "RSB", "resolution": 750},
    "M03": {"wavelength": 0.488, "type": "RSB", "resolution": 750},
    "M04": {"wavelength": 0.555, "type": "RSB", "resolution": 750},
    "M05": {"wavelength": 0.672, "type": "RSB", "resolution": 750},
    "M06": {"wavelength": 0.746, "type": "RSB", "resolution": 750},
    "M07": {"wavelength": 0.865, "type": "RSB", "resolution": 750},
    "M08": {"wavelength": 1.240, "type": "RSB", "resolution": 750},
    "M09": {"wavelength": 1.378, "type": "RSB", "resolution": 750},
    "M10": {"wavelength": 1.610, "type": "RSB", "resolution": 750},
    "M11": {"wavelength": 2.250, "type": "RSB", "resolution": 750},
    "M12": {"wavelength": 3.700, "type": "TEB", "resolution": 750},
    "M13": {"wavelength": 4.050, "type": "TEB", "resolution": 750},
    "M14": {"wavelength": 8.550, "type": "TEB", "resolution": 750},
    "M15": {"wavelength": 10.763, "type": "TEB", "resolution": 750},
    "M16": {"wavelength": 12.013, "type": "TEB", "resolution": 750},
}

VNP21_LAYERS = {
    "LST": {"description": "Land Surface Temperature", "units": "Kelvin", "wavelength": None},
    "Emis_14": {"description": "Emissivity Band M14", "units": "dimensionless", "wavelength": 8.550},
    "Emis_15": {"description": "Emissivity Band M15", "units": "dimensionless", "wavelength": 10.763},
    "Emis_16": {"description": "Emissivity Band M16", "units": "dimensionless", "wavelength": 12.013},
    "LST_err": {"description": "LST Error", "units": "Kelvin", "wavelength": None},
    "Emis_14_err": {"description": "Emissivity M14 Error", "units": "dimensionless", "wavelength": None},
    "Emis_15_err": {"description": "Emissivity M15 Error", "units": "dimensionless", "wavelength": None},
    "Emis_16_err": {"description": "Emissivity M16 Error", "units": "dimensionless", "wavelength": None},
    "QC": {"description": "Quality Control", "units": "bit field", "wavelength": None},
    "View_angle": {"description": "View Zenith Angle", "units": "degrees", "wavelength": None},
    "Emis_ASTER": {"description": "ASTER GED Emissivity", "units": "dimensionless", "wavelength": None},
    "PWV": {"description": "Precipitable Water Vapor", "units": "cm", "wavelength": None},
}

# Filename patterns for product identification
PRODUCT_PATTERNS = {
    "VNP02IMG": re.compile(r"V[NJ][P1]02IMG", re.IGNORECASE),
    "VNP02MOD": re.compile(r"V[NJ][P1]02MOD", re.IGNORECASE),
    "VNP21": re.compile(r"V[NJ][P1]21\.", re.IGNORECASE),
}

GEO_PATTERNS = {
    "VNP02IMG": re.compile(r"V[NJ][P1]03IMG", re.IGNORECASE),
    "VNP02MOD": re.compile(r"V[NJ][P1]03MOD", re.IGNORECASE),
}


# =============================================================================
# Utility Functions
# =============================================================================

def identify_product(filepath):
    """Identify the VIIRS product type from the filename."""
    basename = os.path.basename(filepath)
    for ptype, pattern in PRODUCT_PATTERNS.items():
        if pattern.search(basename):
            return ptype
    return None


def find_geo_file(data_filepath, product_type):
    """Auto-detect matching VNP03 geolocation file in the same directory."""
    if product_type == "VNP21":
        return None  # VNP21 has embedded geolocation

    directory = os.path.dirname(data_filepath)
    data_basename = os.path.basename(data_filepath)

    # Extract the acquisition timestamp: AYYYYDDD.HHMM
    match = re.search(r"\.A(\d{7})\.(\d{4})\.", data_basename)
    if not match:
        return None

    date_str = match.group(1)
    time_str = match.group(2)
    geo_pattern = GEO_PATTERNS.get(product_type)
    if not geo_pattern:
        return None

    # Search for matching geo file in the same directory
    for f in os.listdir(directory):
        if geo_pattern.search(f):
            geo_match = re.search(r"\.A(\d{7})\.(\d{4})\.", f)
            if geo_match and geo_match.group(1) == date_str and geo_match.group(2) == time_str:
                return os.path.join(directory, f)
    return None


def get_utm_epsg(center_lon, center_lat):
    """Determine the UTM zone EPSG code from center coordinates."""
    zone = int((center_lon + 180) / 6) + 1
    if center_lat >= 0:
        return 32600 + zone
    else:
        return 32700 + zone


def read_viirs_data(filepath, product_type, selected_bands, geo_filepath=None, log_fn=None):
    """
    Read VIIRS data and geolocation from netCDF files.
    All L1B bands are output as scaled radiance (W/m^2/um/sr).
    Returns: dict of {band_name: 2D_array}, lat_array, lon_array, metadata_dict
    """
    if log_fn is None:
        log_fn = print

    ds = nc4.Dataset(filepath, "r")
    data = {}
    metadata = {"product_type": product_type, "source_file": os.path.basename(filepath)}

    # Read global attributes
    for attr in ["platform", "instrument", "StartTime", "EndTime", "DayNightFlag"]:
        try:
            metadata[attr] = ds.getncattr(attr)
        except AttributeError:
            pass

    if product_type in ("VNP02IMG", "VNP02MOD"):
        # L1B products: data in observation_data group
        obs_group = ds.groups.get("observation_data")
        if obs_group is None:
            obs_group = ds

        for band_name in selected_bands:
            log_fn(f"  Reading band {band_name}...")
            if band_name not in obs_group.variables:
                log_fn(f"  WARNING: Band {band_name} not found in file, skipping.")
                continue

            var = obs_group.variables[band_name]

            # Disable auto-scaling so we get raw unsigned integers
            var.set_auto_maskandscale(False)
            raw = var[:]

            # Get radiance scale factors and validity range
            rad_scale = getattr(var, "radiance_scale_factor", None)
            rad_offset = getattr(var, "radiance_add_offset", 0.0)
            valid_max = getattr(var, "valid_max", 65527)

            # Build validity mask
            valid = (raw <= valid_max)

            # Apply radiance scaling: radiance = raw * scale + offset
            if rad_scale is not None:
                result = np.full(raw.shape, np.nan, dtype=np.float32)
                result[valid] = raw[valid].astype(np.float32) * float(rad_scale) + float(rad_offset)
                data[band_name] = result
            else:
                # Fallback to generic scale_factor / add_offset
                sf = getattr(var, "scale_factor", 1.0)
                ao = getattr(var, "add_offset", 0.0)
                result = np.full(raw.shape, np.nan, dtype=np.float32)
                result[valid] = raw[valid].astype(np.float32) * float(sf) + float(ao)
                data[band_name] = result

            log_fn(f"  {band_name} radiance range: {np.nanmin(result):.6f} - {np.nanmax(result):.6f} W/m^2/um/sr")

        # Read geolocation from VNP03 file
        if geo_filepath:
            log_fn(f"  Reading geolocation from {os.path.basename(geo_filepath)}...")
            geo_ds = nc4.Dataset(geo_filepath, "r")
            geo_group = geo_ds.groups.get("geolocation_data", geo_ds)
            lat = geo_group.variables["latitude"][:]
            lon = geo_group.variables["longitude"][:]
            geo_ds.close()
        else:
            raise ValueError("No geolocation file found for L1B product.")

    elif product_type == "VNP21":
        # VNP21: data may be at root level or inside a group (e.g., VIIRS_Swath_LSTE)
        # Find the correct group containing the science data
        data_source = ds
        if len(ds.variables) == 0 and len(ds.groups) > 0:
            # Data is inside a group — search for it
            for grp_name in ds.groups:
                grp = ds.groups[grp_name]
                # Check if this group has our variables directly
                if "LST" in grp.variables:
                    data_source = grp
                    log_fn(f"  Data found in group: {grp_name}")
                    break
                # Check subgroups (HDF-EOS5 often nests Data_Fields and Geolocation_Fields)
                for sub_name in grp.groups:
                    sub = grp.groups[sub_name]
                    if "LST" in sub.variables:
                        data_source = sub
                        log_fn(f"  Data found in group: {grp_name}/{sub_name}")
                        break
                if data_source is not ds:
                    break

        # Read science data layers
        for layer_name in selected_bands:
            log_fn(f"  Reading layer {layer_name}...")

            # Search for the variable in data_source and its subgroups
            var = None
            if layer_name in data_source.variables:
                var = data_source.variables[layer_name]
            else:
                # Search all subgroups
                for sg_name in data_source.groups:
                    sg = data_source.groups[sg_name]
                    if layer_name in sg.variables:
                        var = sg.variables[layer_name]
                        break

            if var is None:
                log_fn(f"  WARNING: Layer {layer_name} not found, skipping.")
                continue

            raw = var[:]

            # Convert masked array to float32 with NaN for masked values
            if hasattr(raw, "filled"):
                result = raw.filled(np.nan).astype(np.float32)
            else:
                result = raw.astype(np.float32)

            data[layer_name] = result
            valid_count = np.sum(np.isfinite(result))
            log_fn(f"  {layer_name}: {valid_count:,} valid pixels, "
                   f"range: {np.nanmin(result):.4f} - {np.nanmax(result):.4f}")

        # Find geolocation — search everywhere in the file
        lat_var = None
        lon_var = None

        def find_var_recursive(group, names):
            """Search for a variable by name in a group and all its subgroups."""
            for name in names:
                if name in group.variables:
                    return group.variables[name]
            for sg_name in group.groups:
                result = find_var_recursive(group.groups[sg_name], names)
                if result is not None:
                    return result
            return None

        lat_var = find_var_recursive(ds, ["Latitude", "latitude", "lat"])
        lon_var = find_var_recursive(ds, ["Longitude", "longitude", "lon"])

        if lat_var is None or lon_var is None:
            ds.close()
            # Build a full listing of all variables for debugging
            def list_all_vars(group, prefix=""):
                items = []
                for v in group.variables:
                    items.append(f"{prefix}{v}")
                for g in group.groups:
                    items.extend(list_all_vars(group.groups[g], f"{prefix}{g}/"))
                return items
            all_vars = list_all_vars(ds)
            raise ValueError(f"Could not find latitude/longitude in VNP21 file. "
                           f"All variables: {all_vars}")

        lat = lat_var[:]
        lon = lon_var[:]

        # Convert masked arrays
        if hasattr(lat, "filled"):
            lat = lat.filled(np.nan)
        if hasattr(lon, "filled"):
            lon = lon.filled(np.nan)

        log_fn(f"  Geolocation: lat [{np.nanmin(lat):.2f} to {np.nanmax(lat):.2f}], "
               f"lon [{np.nanmin(lon):.2f} to {np.nanmax(lon):.2f}]")
    else:
        ds.close()
        raise ValueError(f"Unsupported product type: {product_type}")

    ds.close()

    # Convert masked arrays to regular arrays with NaN
    if hasattr(lat, "filled"):
        lat = lat.filled(np.nan)
    if hasattr(lon, "filled"):
        lon = lon.filled(np.nan)
    for k in data:
        if hasattr(data[k], "filled"):
            data[k] = data[k].filled(np.nan)

    return data, lat.astype(np.float64), lon.astype(np.float64), metadata


def reproject_swath(data_dict, lat, lon, projection, resampling, pixel_size_m, log_fn=None):
    """
    Reproject swath data to a regular grid using KD-tree for fast lookup.
    Returns: dict of {band: 2D_gridded_array}, transform, crs, (rows, cols)
    """
    if log_fn is None:
        log_fn = print

    import time
    t_start = time.time()

    # Build a common validity mask across geolocation
    valid_mask = np.isfinite(lat) & np.isfinite(lon)
    n_valid = np.sum(valid_mask)
    log_fn(f"  Valid source pixels: {n_valid:,} of {lat.size:,}")

    if n_valid < 10:
        log_fn("  ERROR: Too few valid geolocation pixels.")
        return {}, None, None, (0, 0)

    center_lat = float(np.nanmean(lat[valid_mask]))
    center_lon = float(np.nanmean(lon[valid_mask]))

    if projection == "UTM":
        epsg = get_utm_epsg(center_lon, center_lat)
        out_crs = CRS.from_epsg(epsg)
        log_fn(f"  UTM Zone: EPSG:{epsg}")

        transformer = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)
        log_fn("  Transforming coordinates to UTM...")
        x_all, y_all = transformer.transform(lon[valid_mask], lat[valid_mask])
        pixel_size = pixel_size_m
    else:
        out_crs = CRS.from_epsg(4326)
        x_all = lon[valid_mask].astype(np.float64)
        y_all = lat[valid_mask].astype(np.float64)
        pixel_size = pixel_size_m / 111320.0

    # Determine grid bounds
    x_min, x_max = float(np.min(x_all)), float(np.max(x_all))
    y_min, y_max = float(np.min(y_all)), float(np.max(y_all))

    # Add small buffer
    x_min -= pixel_size
    x_max += pixel_size
    y_min -= pixel_size
    y_max += pixel_size

    # Create output grid
    cols = int(np.ceil((x_max - x_min) / pixel_size))
    rows = int(np.ceil((y_max - y_min) / pixel_size))

    # Cap grid size to prevent memory issues
    max_dim = 15000
    if cols > max_dim or rows > max_dim:
        scale = max(cols, rows) / max_dim
        pixel_size *= scale
        cols = int(np.ceil((x_max - x_min) / pixel_size))
        rows = int(np.ceil((y_max - y_min) / pixel_size))
        log_fn(f"  Grid capped to {cols}x{rows} (pixel size adjusted to {pixel_size:.6f})")

    log_fn(f"  Output grid: {cols} x {rows} pixels ({cols * rows:,} total)")

    # Subsample source points if very large (>2M points)
    # Use strided subsampling that preserves edge pixels
    max_source_pts = 2_000_000
    if n_valid > max_source_pts:
        step = int(np.ceil(n_valid / max_source_pts))
        log_fn(f"  Subsampling source points: every {step}th pixel ({n_valid:,} -> ~{n_valid // step:,})")

        # Get valid indices in the original flattened array
        all_valid_indices = np.where(valid_mask.ravel())[0]

        # Subsample the middle but keep all edge rows/columns
        # Edge pixels are critical since they're already sparse
        n_rows, n_cols = lat.shape
        edge_width = max(4, n_cols // 20)  # ~5% from each side

        # Build mask of edge pixels (left/right edges of each scan line)
        row_indices, col_indices = np.unravel_index(all_valid_indices, (n_rows, n_cols))
        is_edge = (col_indices < edge_width) | (col_indices >= n_cols - edge_width)

        # Keep all edge pixels, subsample the rest
        edge_idx = all_valid_indices[is_edge]
        interior_idx = all_valid_indices[~is_edge][::step]
        combined_idx = np.sort(np.concatenate([edge_idx, interior_idx]))

        valid_indices = combined_idx
        x_src = np.empty(len(combined_idx), dtype=np.float64)
        y_src = np.empty(len(combined_idx), dtype=np.float64)

        # Map combined indices back to x_all/y_all positions
        # all_valid_indices maps position in valid set -> flat array index
        # We need the reverse: flat array index -> position in x_all/y_all
        idx_lookup = np.full(lat.size, -1, dtype=np.int64)
        idx_lookup[all_valid_indices] = np.arange(n_valid)
        positions = idx_lookup[combined_idx]

        x_src = x_all[positions]
        y_src = y_all[positions]

        log_fn(f"  Kept {len(edge_idx):,} edge + {len(interior_idx):,} interior = {len(combined_idx):,} source points")
    else:
        x_src = x_all
        y_src = y_all
        valid_indices = np.where(valid_mask.ravel())[0]

    n_src = len(x_src)
    log_fn(f"  Building KD-tree with {n_src:,} source points...")
    t_tree = time.time()

    src_points = np.column_stack((x_src, y_src))
    tree = cKDTree(src_points)

    log_fn(f"  KD-tree built in {time.time() - t_tree:.1f}s")

    # Build output grid coordinates (pixel centers)
    grid_x = np.linspace(x_min + pixel_size / 2, x_max - pixel_size / 2, cols)
    grid_y = np.linspace(y_max - pixel_size / 2, y_min + pixel_size / 2, rows)
    grid_xx, grid_yy = np.meshgrid(grid_x, grid_y)
    grid_pts = np.column_stack((grid_xx.ravel(), grid_yy.ravel()))

    # Query KD-tree: find nearest source pixel for each grid cell
    log_fn("  Querying nearest neighbors...")
    t_query = time.time()

    # Set max distance: larger at edges due to pixel footprint expansion
    # VIIRS pixels at edge of scan can be ~1.6km for 750m nadir bands
    # Use 5x pixel size to accommodate the bowtie distortion at scan edges
    max_dist = pixel_size * 5.0
    distances, indices = tree.query(grid_pts, k=1, workers=-1)

    log_fn(f"  Query complete in {time.time() - t_query:.1f}s")

    # Create distance mask (reject grid cells too far from any source pixel)
    dist_mask = distances <= max_dist
    n_filled = np.sum(dist_mask)
    log_fn(f"  Grid cells with data: {n_filled:,} of {len(grid_pts):,} ({100 * n_filled / len(grid_pts):.1f}%)")

    if resampling == "Nearest Neighbor":
        # Simple nearest neighbor: just grab the value at the nearest source pixel
        out_data = {}
        for band_name, band_arr in data_dict.items():
            log_fn(f"  Mapping {band_name} (nearest neighbor)...")
            flat_band = band_arr.ravel()

            result = np.full(len(grid_pts), np.nan, dtype=np.float32)
            valid_grid = dist_mask
            src_idx = valid_indices[indices[valid_grid]]
            result[valid_grid] = flat_band[src_idx]

            out_data[band_name] = result.reshape(rows, cols)
            valid_count = np.sum(np.isfinite(out_data[band_name]))
            log_fn(f"  {band_name}: {valid_count:,} valid output pixels")
    else:
        # Bilinear: use KD-tree with k=4 neighbors and inverse distance weighting
        log_fn("  Querying 4 nearest neighbors for bilinear interpolation...")
        t_bilin = time.time()
        distances_4, indices_4 = tree.query(grid_pts, k=4, workers=-1)
        log_fn(f"  4-neighbor query in {time.time() - t_bilin:.1f}s")

        out_data = {}
        for band_name, band_arr in data_dict.items():
            log_fn(f"  Mapping {band_name} (IDW bilinear)...")
            flat_band = band_arr.ravel()

            result = np.full(len(grid_pts), np.nan, dtype=np.float32)

            # For each grid point, compute inverse-distance-weighted average
            # Only use grid points where the nearest neighbor is within max_dist
            valid_grid = dist_mask

            if np.sum(valid_grid) > 0:
                d = distances_4[valid_grid]
                idx = indices_4[valid_grid]

                # Get values at neighbor locations
                vals = np.full_like(d, np.nan, dtype=np.float32)
                for k in range(4):
                    src_idx = valid_indices[idx[:, k]]
                    vals[:, k] = flat_band[src_idx]

                # Handle exact matches (distance = 0)
                exact = d < 1e-10
                has_exact = np.any(exact, axis=1)

                # IDW weights
                with np.errstate(divide="ignore", invalid="ignore"):
                    weights = 1.0 / d
                weights[~np.isfinite(vals)] = 0.0
                weights[~np.isfinite(weights)] = 0.0

                w_sum = np.sum(weights, axis=1)
                v_weighted = np.nansum(weights * vals, axis=1)

                with np.errstate(divide="ignore", invalid="ignore"):
                    interpolated = v_weighted / w_sum

                # For exact matches, use the exact value
                for k in range(4):
                    mask = has_exact & exact[:, k]
                    interpolated[mask] = vals[mask, k]

                interpolated[w_sum == 0] = np.nan
                result[valid_grid] = interpolated

            out_data[band_name] = result.reshape(rows, cols)
            valid_count = np.sum(np.isfinite(out_data[band_name]))
            log_fn(f"  {band_name}: {valid_count:,} valid output pixels")

    transform = from_bounds(x_min, y_min, x_max, y_max, cols, rows)

    # --- Gap-fill pass for bowtie deletion regions ---
    # VIIRS deletes redundant pixels at scan edges (bowtie deletion):
    #   0 deleted at <31.59°, 2 deleted at 31.59-44.68°, 4 deleted at >44.68°
    # These show up as gaps after reprojection. Fill them using iterative
    # neighbor averaging to close thin gaps while preserving real edges.
    log_fn("  Filling bowtie gaps...")
    t_fill = time.time()

    for band_name in list(out_data.keys()):
        arr = out_data[band_name]
        total_pixels = arr.size
        valid_before = np.sum(np.isfinite(arr))

        if valid_before == 0 or valid_before == total_pixels:
            continue

        filled = arr.copy()

        # 5 passes: first 3 require 2+ neighbors, last 2 require 1+ neighbor
        # This fills wider gaps progressively
        for fill_pass in range(5):
            current_nan = np.isnan(filled)
            if not np.any(current_nan):
                break

            min_neighbors = 2 if fill_pass < 3 else 1

            # Check 4 cardinal neighbors plus 4 diagonal neighbors
            padded = np.pad(filled, 1, mode="constant", constant_values=np.nan)
            neighbors = np.stack([
                padded[0:-2, 1:-1],  # above
                padded[2:,   1:-1],  # below
                padded[1:-1, 0:-2],  # left
                padded[1:-1, 2:],    # right
                padded[0:-2, 0:-2],  # top-left
                padded[0:-2, 2:],    # top-right
                padded[2:,   0:-2],  # bottom-left
                padded[2:,   2:],    # bottom-right
            ], axis=0)

            neighbor_valid = np.isfinite(neighbors)
            n_valid_neighbors = np.sum(neighbor_valid, axis=0)

            fill_mask = current_nan & (n_valid_neighbors >= min_neighbors)

            if not np.any(fill_mask):
                break

            with np.errstate(invalid="ignore"):
                neighbor_sum = np.nansum(neighbors, axis=0)
                neighbor_mean = neighbor_sum / n_valid_neighbors

            filled[fill_mask] = neighbor_mean[fill_mask]

        valid_after = np.sum(np.isfinite(filled))
        n_filled = valid_after - valid_before
        if n_filled > 0:
            out_data[band_name] = filled
            log_fn(f"  {band_name}: filled {n_filled:,} gap pixels")

    log_fn(f"  Gap-fill complete in {time.time() - t_fill:.1f}s")
    log_fn(f"  Total reprojection time: {time.time() - t_start:.1f}s")
    return out_data, transform, out_crs, (rows, cols)


def get_wavelength(band_name, product_type):
    """Get the wavelength in micrometers for a given band/layer name."""
    if product_type in ("VNP02IMG", "VNP02MOD"):
        info = VIIRS_BAND_INFO.get(band_name, {})
        return info.get("wavelength", 0.0)
    else:
        info = VNP21_LAYERS.get(band_name, {})
        return info.get("wavelength", 0.0) or 0.0


def write_geotiff(out_path, data_dict, transform, crs, shape, product_type,
                  band_order=None, log_fn=None):
    """Write reprojected data to a multi-band GeoTIFF with wavelength metadata."""
    if log_fn is None:
        log_fn = print

    if band_order is None:
        band_order = list(data_dict.keys())

    n_bands = len(band_order)
    rows, cols = shape

    profile = {
        "driver": "GTiff",
        "dtype": "float32",
        "width": cols,
        "height": rows,
        "count": n_bands,
        "crs": crs,
        "transform": transform,
        "nodata": np.nan,
        "compress": "lzw",
        "tiled": True,
        "blockxsize": 256,
        "blockysize": 256,
    }

    with rasterio.open(out_path, "w", **profile) as dst:
        for i, band_name in enumerate(band_order, 1):
            dst.write(data_dict[band_name], i)
            wl = get_wavelength(band_name, product_type)
            desc = f"{band_name} ({wl:.3f} um)" if wl > 0 else band_name
            dst.set_band_description(i, desc)

        # Write wavelengths as dataset-level metadata
        wavelength_str = ", ".join(
            f"{get_wavelength(b, product_type):.3f}" for b in band_order
        )
        dst.update_tags(
            WAVELENGTHS=wavelength_str,
            WAVELENGTH_UNITS="Micrometers",
            BAND_NAMES=", ".join(band_order),
            PRODUCT_TYPE=product_type,
        )

    log_fn(f"  GeoTIFF written: {os.path.basename(out_path)}")


def write_envi(out_path, data_dict, transform, crs, shape, product_type,
               band_order=None, log_fn=None):
    """Write reprojected data to ENVI format with full control over the header."""
    if log_fn is None:
        log_fn = print

    if band_order is None:
        band_order = list(data_dict.keys())

    n_bands = len(band_order)
    rows, cols = shape

    # Ensure clean paths
    if out_path.lower().endswith(".hdr"):
        out_path = out_path[:-4]

    dat_path = out_path
    hdr_path = out_path + ".hdr"

    # Write binary data file (BSQ interleave, float32)
    with open(dat_path, "wb") as f:
        for band_name in band_order:
            arr = data_dict[band_name].astype(np.float32)
            f.write(arr.tobytes())

    # Build map info string from the transform
    # Format: {projection, ref_x, ref_y, easting, northing, x_size, y_size, ...}
    x_origin = transform.c  # upper-left x
    y_origin = transform.f  # upper-left y
    x_size = abs(transform.a)
    y_size = abs(transform.e)

    crs_str = crs.to_string() if crs else ""

    if crs and crs.is_geographic:
        map_info = (f"{{Geographic Lat/Lon, 1, 1, {x_origin}, {y_origin}, "
                    f"{x_size}, {y_size}, WGS-84}}")
        coord_sys = ('{GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",'
                     'SPHEROID["WGS_1984",6378137.0,298.257223563]],'
                     'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]}')
    elif crs:
        # UTM
        epsg = crs.to_epsg()
        zone = epsg % 100 if epsg else 0
        hemi = "North" if epsg and epsg < 32700 else "South"
        map_info = (f"{{UTM, 1, 1, {x_origin}, {y_origin}, "
                    f"{x_size}, {y_size}, {zone}, {hemi}, WGS-84}}")
        coord_sys = crs.to_wkt()
        coord_sys = "{" + coord_sys + "}"
    else:
        map_info = None
        coord_sys = None

    # Build wavelength and band name strings
    wavelengths = [f"{get_wavelength(bn, product_type):.4f}" for bn in band_order]
    has_wavelengths = any(float(w) > 0 for w in wavelengths)

    # Determine byte order (system native)
    import sys as _sys
    byte_order = 0 if _sys.byteorder == "little" else 1

    # Write the header file from scratch
    with open(hdr_path, "w") as f:
        f.write("ENVI\n")
        f.write(f"description = {{{os.path.basename(out_path)}}}\n")
        f.write(f"samples = {cols}\n")
        f.write(f"lines = {rows}\n")
        f.write(f"bands = {n_bands}\n")
        f.write("header offset = 0\n")
        f.write("file type = ENVI Standard\n")
        f.write("data type = 4\n")
        f.write("interleave = bsq\n")
        f.write(f"byte order = {byte_order}\n")

        if map_info:
            f.write(f"map info = {map_info}\n")
        if coord_sys:
            f.write(f"coordinate system string = {coord_sys}\n")

        # Wavelengths
        if has_wavelengths:
            wl_str = ", ".join(wavelengths)
            f.write(f"wavelength = {{{wl_str}}}\n")
            f.write("wavelength units = Micrometers\n")

        # Band names
        bn_str = ", ".join(band_order)
        f.write(f"band names = {{{bn_str}}}\n")

        # NoData
        f.write("data ignore value = nan\n")

    log_fn(f"  ENVI written: {os.path.basename(out_path)}")

    log_fn(f"  ENVI written: {os.path.basename(out_path)}")


# =============================================================================
# GUI Application
# =============================================================================

class VIIRSConverterApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("VIIRS Data Converter v1.0")
        self.geometry("920x950")
        self.minsize(880, 900)

        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")

        self.input_files = []
        self.processing = False

        self._build_ui()
        self._check_dependencies()

    def _build_ui(self):
        # Main container
        main = ctk.CTkFrame(self, fg_color="transparent")
        main.pack(fill="both", expand=True, padx=10, pady=10)
        main.grid_columnconfigure(0, weight=1)

        # Title
        title_frame = ctk.CTkFrame(main, fg_color="#1a1a2e", corner_radius=10)
        title_frame.grid(row=0, column=0, sticky="ew", pady=(0, 10))
        ctk.CTkLabel(
            title_frame,
            text="VIIRS Data Converter v1.0",
            font=ctk.CTkFont(size=22, weight="bold"),
            text_color="#4fc3f7",
        ).pack(pady=(10, 2))
        ctk.CTkLabel(
            title_frame,
            text="Convert VIIRS swath products to GeoTIFF & ENVI",
            font=ctk.CTkFont(size=12),
            text_color="#90a4ae",
        ).pack(pady=(0, 10))

        # --- Input Section ---
        input_frame = ctk.CTkFrame(main, fg_color="#1e1e2e", corner_radius=8)
        input_frame.grid(row=1, column=0, sticky="ew", pady=(0, 8))
        input_frame.grid_columnconfigure(0, weight=1)

        ctk.CTkLabel(
            input_frame, text="Input Files", font=ctk.CTkFont(size=14, weight="bold")
        ).grid(row=0, column=0, sticky="w", padx=12, pady=(10, 4))

        btn_row = ctk.CTkFrame(input_frame, fg_color="transparent")
        btn_row.grid(row=1, column=0, sticky="ew", padx=12, pady=(0, 4))

        ctk.CTkButton(
            btn_row, text="Select File(s)", width=130, command=self._select_files
        ).pack(side="left", padx=(0, 8))
        ctk.CTkButton(
            btn_row, text="Select Folder", width=130, command=self._select_folder
        ).pack(side="left", padx=(0, 8))
        ctk.CTkButton(
            btn_row, text="Clear", width=80, fg_color="#555", command=self._clear_files
        ).pack(side="left")

        self.file_label = ctk.CTkLabel(
            input_frame, text="No files selected", text_color="#90a4ae", anchor="w"
        )
        self.file_label.grid(row=2, column=0, sticky="ew", padx=12, pady=(0, 10))

        # --- Band Selection ---
        band_frame = ctk.CTkFrame(main, fg_color="#1e1e2e", corner_radius=8)
        band_frame.grid(row=2, column=0, sticky="ew", pady=(0, 8))
        band_frame.grid_columnconfigure(0, weight=1)

        ctk.CTkLabel(
            band_frame,
            text="Band / Layer Selection",
            font=ctk.CTkFont(size=14, weight="bold"),
        ).grid(row=0, column=0, sticky="w", padx=12, pady=(10, 4))

        self.band_info_label = ctk.CTkLabel(
            band_frame,
            text="Select input files to see available bands",
            text_color="#90a4ae",
            anchor="w",
        )
        self.band_info_label.grid(row=1, column=0, sticky="ew", padx=12, pady=(0, 4))

        self.band_scroll = ctk.CTkScrollableFrame(
            band_frame, height=60, fg_color="#16162a"
        )
        self.band_scroll.grid(row=2, column=0, sticky="ew", padx=12, pady=(0, 4))
        self.band_scroll.grid_columnconfigure((0, 1, 2, 3, 4), weight=1)

        band_btn_row = ctk.CTkFrame(band_frame, fg_color="transparent")
        band_btn_row.grid(row=3, column=0, sticky="w", padx=12, pady=(0, 10))
        ctk.CTkButton(
            band_btn_row, text="Select All", width=90, command=self._select_all_bands
        ).pack(side="left", padx=(0, 8))
        ctk.CTkButton(
            band_btn_row, text="Deselect All", width=90, fg_color="#555", command=self._deselect_all_bands
        ).pack(side="left")

        self.band_vars = {}

        # --- Output Options (compact two-column layout) ---
        opts_frame = ctk.CTkFrame(main, fg_color="#1e1e2e", corner_radius=8)
        opts_frame.grid(row=3, column=0, sticky="ew", pady=(0, 8))

        ctk.CTkLabel(
            opts_frame,
            text="Output Options",
            font=ctk.CTkFont(size=14, weight="bold"),
        ).pack(anchor="w", padx=12, pady=(10, 6))

        cols_row = ctk.CTkFrame(opts_frame, fg_color="transparent")
        cols_row.pack(fill="x", padx=12, pady=(0, 10))

        # Left column: Projection, Resampling, Format
        left_col = ctk.CTkFrame(cols_row, fg_color="transparent")
        left_col.pack(side="left", fill="both", expand=True, padx=(0, 8))

        ctk.CTkLabel(left_col, text="Projection:", font=ctk.CTkFont(size=12)).pack(anchor="w")
        self.proj_var = ctk.StringVar(value="Geographic (EPSG:4326)")
        ctk.CTkOptionMenu(
            left_col, variable=self.proj_var,
            values=["Geographic (EPSG:4326)", "UTM (Auto-detect)"], width=200,
        ).pack(anchor="w", pady=(0, 6))

        ctk.CTkLabel(left_col, text="Resampling:", font=ctk.CTkFont(size=12)).pack(anchor="w")
        self.resample_var = ctk.StringVar(value="Nearest Neighbor")
        ctk.CTkOptionMenu(
            left_col, variable=self.resample_var,
            values=["Nearest Neighbor", "Bilinear"], width=200,
        ).pack(anchor="w", pady=(0, 6))

        ctk.CTkLabel(left_col, text="Output Format:", font=ctk.CTkFont(size=12)).pack(anchor="w")
        self.fmt_var = ctk.StringVar(value="Both (GeoTIFF + ENVI)")
        ctk.CTkOptionMenu(
            left_col, variable=self.fmt_var,
            values=["GeoTIFF Only", "ENVI Only", "Both (GeoTIFF + ENVI)"], width=200,
        ).pack(anchor="w")

        # Right column: Output Dir, Overwrite, Convert button, Progress
        right_col = ctk.CTkFrame(cols_row, fg_color="transparent")
        right_col.pack(side="left", fill="both", expand=True, padx=(8, 0))

        ctk.CTkLabel(right_col, text="Output Directory:", font=ctk.CTkFont(size=12)).pack(anchor="w")
        outdir_row = ctk.CTkFrame(right_col, fg_color="transparent")
        outdir_row.pack(fill="x", pady=(0, 6))
        self.outdir_var = ctk.StringVar(value="Same as input")
        ctk.CTkEntry(outdir_row, textvariable=self.outdir_var, width=160).pack(
            side="left", fill="x", expand=True, padx=(0, 4)
        )
        ctk.CTkButton(
            outdir_row, text="Browse", width=70, command=self._browse_outdir
        ).pack(side="left")

        self.overwrite_var = ctk.BooleanVar(value=False)
        ctk.CTkCheckBox(
            right_col, text="Overwrite existing files", variable=self.overwrite_var,
        ).pack(anchor="w", pady=(4, 6))

        self.process_btn = ctk.CTkButton(
            right_col,
            text="Convert",
            font=ctk.CTkFont(size=13, weight="bold"),
            height=32,
            command=self._start_processing,
        )
        self.process_btn.pack(fill="x", pady=(2, 0))

        self.progress = ctk.CTkProgressBar(right_col, height=6)
        self.progress.pack(fill="x", pady=(6, 0))
        self.progress.set(0)

        # --- Log ---
        log_frame = ctk.CTkFrame(main, fg_color="#1e1e2e", corner_radius=8)
        log_frame.grid(row=4, column=0, sticky="nsew", pady=(0, 0))
        log_frame.grid_columnconfigure(0, weight=1)
        log_frame.grid_rowconfigure(1, weight=1)
        main.grid_rowconfigure(4, weight=1)

        ctk.CTkLabel(
            log_frame, text="Processing Log", font=ctk.CTkFont(size=14, weight="bold")
        ).grid(row=0, column=0, sticky="w", padx=12, pady=(10, 4))

        self.log_text = ctk.CTkTextbox(
            log_frame, height=100, font=ctk.CTkFont(family="Consolas", size=11),
            fg_color="#0d0d1a", text_color="#c0c0c0"
        )
        self.log_text.grid(row=1, column=0, sticky="nsew", padx=12, pady=(0, 10))

    def _check_dependencies(self):
        """Check for required libraries."""
        missing = []
        if nc4 is None:
            missing.append("netCDF4")
        if rasterio is None:
            missing.append("rasterio")
        if Transformer is None:
            missing.append("pyproj")
        if griddata is None:
            missing.append("scipy")

        if missing:
            self.log(f"WARNING: Missing dependencies: {', '.join(missing)}")
            self.log("Install with: pip install " + " ".join(missing))
            self.process_btn.configure(state="disabled")
        else:
            self.log("All dependencies OK. Ready to convert.")

    def log(self, message):
        """Add a message to the processing log (thread-safe via after)."""
        self.after(0, self._log_to_widget, message)

    def _log_to_widget(self, message):
        """Actually write to the log widget (must run on main thread)."""
        self.log_text.insert("end", message + "\n")
        self.log_text.see("end")

    def _select_files(self):
        files = filedialog.askopenfilenames(
            title="Select VIIRS Data Files",
            filetypes=[
                ("VIIRS Files", "*.nc *.h5 *.hdf5 *.he5"),
                ("NetCDF", "*.nc"),
                ("HDF5", "*.h5 *.hdf5 *.he5"),
                ("All Files", "*.*"),
            ],
        )
        if files:
            self.input_files = list(files)
            self._update_file_display()
            self._update_band_selection()

    def _select_folder(self):
        folder = filedialog.askdirectory(title="Select Folder with VIIRS Files")
        if folder:
            files = []
            for ext in ("*.nc", "*.h5", "*.hdf5", "*.he5"):
                files.extend(glob.glob(os.path.join(folder, ext)))
            data_files = []
            for f in files:
                ptype = identify_product(f)
                if ptype is not None:
                    data_files.append(f)
            self.input_files = sorted(data_files)
            self._update_file_display()
            self._update_band_selection()

    def _clear_files(self):
        self.input_files = []
        self.file_label.configure(text="No files selected")
        self._clear_band_checkboxes()
        self.band_info_label.configure(text="Select input files to see available bands")

    def _update_file_display(self):
        n = len(self.input_files)
        if n == 0:
            self.file_label.configure(text="No files selected")
        elif n == 1:
            self.file_label.configure(text=os.path.basename(self.input_files[0]))
        else:
            self.file_label.configure(text=f"{n} files selected")

    def _clear_band_checkboxes(self):
        for widget in self.band_scroll.winfo_children():
            widget.destroy()
        self.band_vars = {}

    def _update_band_selection(self):
        """Update the band selection panel based on the detected product type."""
        self._clear_band_checkboxes()

        if not self.input_files:
            return

        ptype = identify_product(self.input_files[0])
        if ptype is None:
            self.band_info_label.configure(
                text="Could not identify product type from filename"
            )
            return

        geo_file = find_geo_file(self.input_files[0], ptype)
        geo_status = ""
        if ptype in ("VNP02IMG", "VNP02MOD"):
            if geo_file:
                geo_status = f"  |  Geolocation: {os.path.basename(geo_file)}"
            else:
                geo_status = "  |  WARNING: No matching geolocation file found!"

        self.band_info_label.configure(text=f"Product: {ptype}{geo_status}")

        if ptype == "VNP02IMG":
            bands = ["I01", "I02", "I03", "I04", "I05"]
        elif ptype == "VNP02MOD":
            bands = [f"M{i:02d}" for i in range(1, 17)]
        elif ptype == "VNP21":
            bands = list(VNP21_LAYERS.keys())
        else:
            return

        cols = 5 if ptype != "VNP21" else 4
        for i, band in enumerate(bands):
            var = ctk.BooleanVar(value=True)
            self.band_vars[band] = var

            if ptype in ("VNP02IMG", "VNP02MOD"):
                info = VIIRS_BAND_INFO.get(band, {})
                wl = info.get("wavelength", "?")
                btype = info.get("type", "?")
                label = f"{band} ({wl}um, {btype})"
            else:
                label = f"{band}"

            cb = ctk.CTkCheckBox(
                self.band_scroll,
                text=label,
                variable=var,
                font=ctk.CTkFont(size=11),
                checkbox_width=18,
                checkbox_height=18,
            )
            cb.grid(row=i // cols, column=i % cols, sticky="w", padx=4, pady=2)

    def _select_all_bands(self):
        for var in self.band_vars.values():
            var.set(True)

    def _deselect_all_bands(self):
        for var in self.band_vars.values():
            var.set(False)

    def _browse_outdir(self):
        d = filedialog.askdirectory(title="Select Output Directory")
        if d:
            self.outdir_var.set(d)

    def _start_processing(self):
        if self.processing:
            return
        if not self.input_files:
            messagebox.showwarning("No Input", "Please select input files first.")
            return

        selected = [b for b, v in self.band_vars.items() if v.get()]
        if not selected:
            messagebox.showwarning("No Bands", "Please select at least one band/layer.")
            return

        self.processing = True
        self.process_btn.configure(state="disabled", text="Processing...")
        self.log_text.delete("1.0", "end")
        self.progress.set(0)

        thread = threading.Thread(
            target=self._process_files, args=(selected,), daemon=True
        )
        thread.start()

    def _process_files(self, selected_bands):
        """Main processing loop (runs in background thread)."""
        try:
            total = len(self.input_files)
            success = 0
            errors = 0

            projection = "UTM" if "UTM" in self.proj_var.get() else "Geographic"
            resampling = self.resample_var.get()
            fmt = self.fmt_var.get()
            overwrite = self.overwrite_var.get()
            outdir_setting = self.outdir_var.get()

            self.log(f"Processing {total} file(s)...")
            self.log(f"Projection: {projection} | Resampling: {resampling}")
            self.log(f"Format: {fmt} | Overwrite: {overwrite}")
            self.log(f"Selected bands: {', '.join(selected_bands)}")
            self.log("-" * 60)

            for idx, filepath in enumerate(self.input_files):
                basename = os.path.basename(filepath)
                self.log(f"\n[{idx+1}/{total}] {basename}")

                try:
                    ptype = identify_product(filepath)
                    if ptype is None:
                        self.log(f"  ERROR: Cannot identify product type, skipping.")
                        errors += 1
                        continue

                    if ptype == "VNP02IMG":
                        native_res = 375
                    else:
                        native_res = 750

                    geo_file = None
                    if ptype in ("VNP02IMG", "VNP02MOD"):
                        geo_file = find_geo_file(filepath, ptype)
                        if geo_file is None:
                            self.log(f"  ERROR: No matching geolocation file found, skipping.")
                            errors += 1
                            continue
                        self.log(f"  Geolocation: {os.path.basename(geo_file)}")

                    if ptype == "VNP02IMG":
                        avail = [f"I{i:02d}" for i in range(1, 6)]
                    elif ptype == "VNP02MOD":
                        avail = [f"M{i:02d}" for i in range(1, 17)]
                    else:
                        avail = list(VNP21_LAYERS.keys())

                    bands_to_process = [b for b in selected_bands if b in avail]
                    self.log(f"  Product type: {ptype}")
                    self.log(f"  Bands to process: {', '.join(bands_to_process) if bands_to_process else 'NONE'}")
                    if not bands_to_process:
                        self.log(f"  No selected bands available for this product, skipping.")
                        continue

                    self.log(f"  Reading data ({len(bands_to_process)} bands)...")
                    data, lat, lon, meta = read_viirs_data(
                        filepath, ptype, bands_to_process, geo_file, log_fn=self.log
                    )

                    if not data:
                        self.log(f"  ERROR: No valid data read, skipping.")
                        errors += 1
                        continue

                    self.log(f"  Reprojecting to {projection}...")
                    gridded, transform, crs, shape = reproject_swath(
                        data, lat, lon, projection, resampling, native_res, log_fn=self.log
                    )

                    if not gridded:
                        self.log(f"  ERROR: Reprojection produced no output, skipping.")
                        errors += 1
                        continue

                    if outdir_setting == "Same as input":
                        out_dir = os.path.dirname(filepath)
                    else:
                        out_dir = outdir_setting

                    base_name = os.path.splitext(basename)[0]
                    band_order = [b for b in bands_to_process if b in gridded]

                    if fmt in ("GeoTIFF Only", "Both (GeoTIFF + ENVI)"):
                        tiff_dir = os.path.join(out_dir, "GeoTIFF")
                        os.makedirs(tiff_dir, exist_ok=True)
                        tiff_path = os.path.join(tiff_dir, f"{base_name}.tif")

                        if os.path.exists(tiff_path) and not overwrite:
                            self.log(f"  SKIP (exists): {os.path.basename(tiff_path)}")
                        else:
                            write_geotiff(
                                tiff_path, gridded, transform, crs, shape,
                                ptype, band_order=band_order, log_fn=self.log
                            )

                    if fmt in ("ENVI Only", "Both (GeoTIFF + ENVI)"):
                        envi_dir = os.path.join(out_dir, "ENVI")
                        os.makedirs(envi_dir, exist_ok=True)
                        envi_path = os.path.join(envi_dir, base_name)

                        if os.path.exists(envi_path) and not overwrite:
                            self.log(f"  SKIP (exists): {os.path.basename(envi_path)}")
                        else:
                            write_envi(
                                envi_path, gridded, transform, crs, shape,
                                ptype, band_order=band_order, log_fn=self.log
                            )

                    success += 1

                except Exception as e:
                    self.log(f"  ERROR: {str(e)}")
                    self.log(f"  {traceback.format_exc()}")
                    errors += 1

                self.progress.set((idx + 1) / total)

            self.log("\n" + "=" * 60)
            self.log(f"Complete: {success} succeeded, {errors} failed out of {total}")

        except Exception as e:
            self.log(f"\nFATAL ERROR: {str(e)}")
            self.log(traceback.format_exc())

        finally:
            self.processing = False
            self.process_btn.configure(state="normal", text="Convert")


# =============================================================================
# Standalone process_directory for RSDTK integration
# =============================================================================

# Aliases for write functions (avoid name collision with the write_geotiff/write_envi
# parameter names in process_directory)
write_geotiff_func = write_geotiff
write_envi_func = write_envi


def process_directory(input_dir, output_dir=None,
                      write_geotiff=True, write_envi=True,
                      overwrite=True, selected_bands=None):
    """
    Batch process all VIIRS files in a directory.
    Standalone function for integration with RSDTK (no GUI).

    Parameters
    ----------
    selected_bands : list or None
        List of band names to process (e.g. ['M14', 'M15', 'M16']).
        If None, processes all bands for each product type.
        For VNP02IMG: I01-I05
        For VNP02MOD: M01-M16
        For VNP21: LST, Emis_14, Emis_15, Emis_16, etc.
    """
    if output_dir is None:
        output_dir = os.path.join(input_dir, 'converted')

    # Discover all VIIRS HDF5/NetCDF files
    supported_ext = ('.h5', '.nc', '.hdf', '.he5')
    all_files = []
    for f in os.listdir(input_dir):
        if f.lower().endswith(supported_ext):
            fpath = os.path.join(input_dir, f)
            ptype = identify_product(fpath)
            if ptype is not None:
                all_files.append((fpath, ptype))

    if not all_files:
        print("No VIIRS files found.")
        return

    print(f"Found {len(all_files)} VIIRS file(s):\n")

    success = 0
    errors = 0
    for idx, (filepath, ptype) in enumerate(all_files):
        basename = os.path.basename(filepath)
        print(f"\n[{idx+1}/{len(all_files)}] {basename}")
        print(f"  Product type: {ptype}")

        try:
            # Determine native resolution
            native_res = 375 if ptype == "VNP02IMG" else 750

            # Find geolocation file
            geo_file = None
            if ptype in ("VNP02IMG", "VNP02MOD"):
                geo_file = find_geo_file(filepath, ptype)
                if geo_file is None:
                    print(f"  ERROR: No matching geolocation file found, skipping.")
                    errors += 1
                    continue
                print(f"  Geolocation: {os.path.basename(geo_file)}")

            # Select bands — use user selection or default to all
            if ptype == "VNP02IMG":
                all_bands = [f"I{i:02d}" for i in range(1, 6)]
            elif ptype == "VNP02MOD":
                all_bands = [f"M{i:02d}" for i in range(1, 17)]
            else:
                all_bands = list(VNP21_LAYERS.keys())

            if selected_bands is not None:
                # Filter to only bands that exist for this product type
                bands = [b for b in selected_bands if b in all_bands]
                if not bands:
                    print(f"  WARNING: None of {selected_bands} match "
                          f"available bands {all_bands}, using all bands.")
                    bands = all_bands
            else:
                bands = all_bands

            print(f"  Processing bands: {', '.join(bands)}")

            # Read data
            data, lat, lon, meta = read_viirs_data(
                filepath, ptype, bands, geo_file, log_fn=print
            )
            if not data:
                print(f"  ERROR: No valid data read, skipping.")
                errors += 1
                continue

            # Reproject
            print(f"  Reprojecting to Geographic...")
            gridded, transform, crs, shape = reproject_swath(
                data, lat, lon, "Geographic", "IDW", native_res, log_fn=print
            )
            if not gridded:
                print(f"  ERROR: Reprojection produced no output, skipping.")
                errors += 1
                continue

            base_name = os.path.splitext(basename)[0]
            band_order = [b for b in bands if b in gridded]

            if write_geotiff:
                tiff_dir = os.path.join(output_dir, "GeoTIFF")
                os.makedirs(tiff_dir, exist_ok=True)
                tiff_path = os.path.join(tiff_dir, f"{base_name}.tif")
                if not overwrite and os.path.exists(tiff_path):
                    print(f'  [skip] {os.path.basename(tiff_path)} already exists')
                else:
                    write_geotiff_func(
                        tiff_path, gridded, transform, crs, shape,
                        ptype, band_order=band_order, log_fn=print
                    )

            if write_envi:
                envi_dir = os.path.join(output_dir, "ENVI")
                os.makedirs(envi_dir, exist_ok=True)
                envi_path = os.path.join(envi_dir, base_name)
                envi_dat_check = envi_path + '.dat'
                if not overwrite and os.path.exists(envi_dat_check):
                    print(f'  [skip] {os.path.basename(envi_dat_check)} already exists')
                else:
                    write_envi_func(
                        envi_path, gridded, transform, crs, shape,
                        ptype, band_order=band_order, log_fn=print
                    )

            success += 1

        except Exception as e:
            print(f"  ERROR: {str(e)}")
            import traceback
            traceback.print_exc()
            errors += 1

    print(f"\nComplete: {success} succeeded, {errors} failed out of {len(all_files)}")


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    app = VIIRSConverterApp()
    app.mainloop()