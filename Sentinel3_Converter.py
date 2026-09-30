# Decompiled with PyLingual (https://pylingual.io)
# Internal filename: 'Sentinel3_Converter.py'
# Bytecode version: 3.11a7e (3495)
# Source timestamp: 1970-01-01 00:00:00 UTC (0)

"""
Sentinel-3 Product Converter
===============================
Reads Sentinel-3 .SEN3 product folders (netCDF files) and outputs
multiband GeoTIFF and/or ENVI files on a regular geographic grid.

Supported products:
  SLSTR L1B RBT (SL_1_RBT):
    -> VNIR/SWIR radiance: S1-S6 (500m nadir A-stripe) -> 6-band file
    -> TIR brightness temperature: S7-S9, F1-F2 (1km nadir) -> 5-band file

  SLSTR L2 LST (SL_2_LST):
    -> Land Surface Temperature (1km) + uncertainty -> 2-band file

  OLCI L1B EFR (OL_1_EFR):
    -> TOA radiance: Oa01-Oa21 (300m) -> 21-band file

  Synergy L2 SYN (SY_2_SYN):
    -> Surface reflectance: Syn_Oa01-Oa21 (OLCI, 300m) -> up to 21 bands
    -> Surface reflectance: Syn_S1N-S6N (SLSTR nadir, 300m) -> up to 6 bands

Auto-detects product type from .SEN3 folder name.
Supports S3A, S3B, S3C, S3D platforms.

Swath data is reprojected to a regular geographic (lat/lon WGS-84) grid
using scipy cKDTree nearest-neighbor interpolation.

Dependencies: numpy, rasterio, h5py (or netCDF4), scipy
"""
import os
import re
import sys
import glob
import numpy as np
try:
    import rasterio
    from rasterio.crs import CRS
    from rasterio.transform import from_bounds
except ImportError:
    print('ERROR: rasterio is required.')
    sys.exit(1)
try:
    from scipy.spatial import cKDTree
except ImportError:
    print('ERROR: scipy is required.')
    sys.exit(1)
_NC_BACKEND = None
try:
    import netCDF4
    _NC_BACKEND = 'netcdf4'
except ImportError:
    try:
        import h5py
        _NC_BACKEND = 'h5py'
    except ImportError:
        print('ERROR: netCDF4 or h5py is required.')
        sys.exit(1)
def nc_open(filepath):
    """Open a netCDF/HDF5 file. Returns a wrapper object."""
    if _NC_BACKEND == 'netcdf4':
        return netCDF4.Dataset(filepath, 'r')
    else:
        return h5py.File(filepath, 'r')
def nc_read_var(ds, varname):
    """Read a variable as numpy array (raw, no auto-scaling)."""
    if _NC_BACKEND == 'netcdf4':
        var = ds.variables[varname]
        var.set_auto_maskandscale(False)
        return np.array(var[:])
    else:
        return np.array(ds[varname])
def nc_get_attr(ds, varname, attrname, default=None):
    """Get an attribute from a variable."""
    try:
        if _NC_BACKEND == 'netcdf4':
            return ds.variables[varname].getncattr(attrname)
        else:
            return ds[varname].attrs[attrname]
    except (KeyError, AttributeError):
        return default
def nc_has_var(ds, varname):
    """Check if variable exists."""
    if _NC_BACKEND == 'netcdf4':
        return varname in ds.variables
    else:
        return varname in ds
def nc_close(ds):
    """Close the dataset."""
    ds.close()
WGS84_WKT = 'GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433]]'
def build_wgs84_crs():
    """Build WGS-84 CRS from WKT."""
    return CRS.from_wkt(WGS84_WKT)
PRODUCT_PATTERNS = {'SL_1_RBT': re.compile('S3[A-D]_SL_1_RBT', re.IGNORECASE), 'SL_2_LST': re.compile('S3[A-D]_SL_2_LST', re.IGNORECASE), 'OL_1_EFR': re.compile('S3[A-D]_OL_1_EFR', re.IGNORECASE), 'SY_2_SYN': re.compile('S3[A-D]_SY_2_SYN', re.IGNORECASE)}
SLSTR_VNIR_SWIR = {'bands': ['S1', 'S2', 'S3', 'S4', 'S5', 'S6'], 'wavelengths': [0.555, 0.66, 0.87, 1.375, 1.61, 2.25], 'names': ['S1 - Green (0.555 um)', 'S2 - Red (0.660 um)', 'S3 - NIR (0.870 um)', 'S4 - Cirrus (1.375 um)', 'S5 - SWIR-1 (1.610 um)', 'S6 - SWIR-2 (2.250 um)'], 'grid_suffix': 'an', 'geo_file': 'geodetic_an.nc', 'data_type': 'radiance'}
SLSTR_TIR = {'bands': ['S7', 'S8', 'S9', 'F2'], 'wavelengths': [3.74, 10.85, 12.0, 10.85], 'names': ['S7 - MIR (3.740 um)', 'S8 - TIR-1 (10.850 um)', 'S9 - TIR-2 (12.000 um)', 'F2 - Fire TIR (10.850 um)'], 'grid_suffix': 'in', 'geo_file': 'geodetic_in.nc', 'data_type': 'BT'}
SLSTR_FIRE = {'bands': ['F1'], 'wavelengths': [3.74], 'names': ['F1 - Fire MIR (3.740 um)'], 'grid_suffix': 'fn', 'geo_file': 'geodetic_fn.nc', 'data_type': 'BT'}
OLCI_BANDS = {'bands': [f'Oa{i:02d}' for i in range(1, 22)], 'wavelengths': [0.4, 0.4125, 0.4425, 0.49, 0.51, 0.56, 0.62, 0.665, 0.67375, 0.68125, 0.70875, 0.75375, 0.76125, 0.764375, 0.7675, 0.77875, 0.865, 0.885, 0.9, 0.94, 1.02], 'names': ['Oa01 - Aerosol (400 nm)', 'Oa02 - Yellow Sub (412.5 nm)', 'Oa03 - Chl Absorption (442.5 nm)', 'Oa04 - Chl/Water (490 nm)', 'Oa05 - Chl/Sediment (510 nm)', 'Oa06 - Chl Reference (560 nm)', 'Oa07 - Sediment (620 nm)', 'Oa08 - Chl 2nd Max (665 nm)', 'Oa09 - Chl Fluor Peak (673.75 nm)', 'Oa10 - Chl Fluor Base (681.25 nm)', 'Oa11 - Chl Fluor Base (708.75 nm)', 'Oa12 - O2 Absorption (753.75 nm)', 'Oa13 - O2 Band (761.25 nm)', 'Oa14 - Atm Corr (764.375 nm)', 'Oa15 - O2A Cloud Top (767.5 nm)', 'Oa16 - Atm/Aerosol (778.75 nm)', 'Oa17 - Atm/Aerosol (865 nm)', 'Oa18 - Water Vapour (885 nm)', 'Oa19 - Water Vapour/Veg (900 nm)', 'Oa20 - Water Vapour (940 nm)', 'Oa21 - Atm Corr (1020 nm)'], 'geo_file': 'geo_coordinates.nc'}
_SYN_OLCI_INDICES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 17, 18, 21]
SYN_OLCI_BANDS = {'bands': [f'Syn_Oa{i:02d}' for i in _SYN_OLCI_INDICES], 'wavelengths': [OLCI_BANDS['wavelengths'][i - 1] for i in _SYN_OLCI_INDICES], 'names': [f"Syn {OLCI_BANDS['names'][i - 1]}" for i in _SYN_OLCI_INDICES], 'geo_file': 'geolocation.nc'}
SYN_SLSTR_BANDS = {'bands': ['Syn_S1N', 'Syn_S2N', 'Syn_S3N', 'Syn_S5N', 'Syn_S6N'], 'wavelengths': [0.555, 0.66, 0.87, 1.61, 2.25], 'names': ['Syn S1N - Green (0.555 um)', 'Syn S2N - Red (0.660 um)', 'Syn S3N - NIR (0.870 um)', 'Syn S5N - SWIR-1 (1.610 um)', 'Syn S6N - SWIR-2 (2.250 um)'], 'geo_file': 'geolocation.nc'}
def reproject_swath_to_geographic(lat, lon, data_bands, target_resolution=None):
    """
    Reproject swath data (with per-pixel lat/lon) to a regular geographic grid
    using cKDTree nearest-neighbor interpolation.

    Parameters
    ----------
    lat : 2D array of latitudes
    lon : 2D array of longitudes
    data_bands : list of 2D arrays (same shape as lat/lon)
    target_resolution : float, degrees. If None, estimated from data.

    Returns
    -------
    out_bands : list of 2D arrays on the regular grid
    transform : rasterio Affine transform
    out_height, out_width : grid dimensions
    """
    valid = np.isfinite(lat) & np.isfinite(lon)
    valid &= (lat >= (-90)) & (lat <= 90) & (lon >= (-180)) & (lon <= 180)
    lat_valid = lat[valid]
    lon_valid = lon[valid]
    if len(lat_valid) == 0:
        print('    WARNING: No valid coordinates found.')
        return (None, None, 0, 0)
    else:
        lat_min, lat_max = (float(np.nanmin(lat_valid)), float(np.nanmax(lat_valid)))
        lon_min, lon_max = (float(np.nanmin(lon_valid)), float(np.nanmax(lon_valid)))
        if target_resolution is None:
            nrows, ncols = lat.shape
            mid_row = nrows // 2
            mid_col = ncols // 2
            spacings = []
            for col in [mid_col, ncols // 4, 3 * ncols // 4]:
                if col >= ncols:
                    continue
                else:
                    lat_col = lat[:, col]
                    lon_col = lon[:, col]
                    m = np.isfinite(lat_col) & np.isfinite(lon_col)
                    if m.sum() > 10:
                        dlat = np.diff(lat_col[m])
                        dlon = np.diff(lon_col[m])
                        dist = np.sqrt(dlat ** 2 + dlon ** 2)
                        dist = dist[dist > 1e-06]
                        if len(dist) > 0:
                            spacings.append(float(np.median(dist)))
            for row in [mid_row, nrows // 4, 3 * nrows // 4]:
                if row >= nrows:
                    continue
                else:
                    lat_row = lat[row, :]
                    lon_row = lon[row, :]
                    m = np.isfinite(lat_row) & np.isfinite(lon_row)
                    if m.sum() > 10:
                        dlat = np.diff(lat_row[m])
                        dlon = np.diff(lon_row[m])
                        dist = np.sqrt(dlat ** 2 + dlon ** 2)
                        dist = dist[dist > 1e-06]
                        if len(dist) > 0:
                            spacings.append(float(np.median(dist)))
            if spacings:
                pixel_spacing = float(np.median(spacings))
            else:
                pixel_spacing = 0.005
            lat_extent = lat_max - lat_min
            lon_extent = lon_max - lon_min
            n_total = np.sum(valid)
            if n_total > 0:
                extent_res = np.sqrt(lat_extent * lon_extent / n_total)
            else:
                extent_res = 0.005
            target_resolution = max(pixel_spacing, extent_res)
            target_resolution = max(target_resolution, 1e-05)
        out_height = max(1, int(np.ceil((lat_max - lat_min) / target_resolution)))
        out_width = max(1, int(np.ceil((lon_max - lon_min) / target_resolution)))
        max_dim = 20000
        if out_height > max_dim or out_width > max_dim:
            scale_factor = max(out_height, out_width) / max_dim
            target_resolution *= scale_factor
            out_height = max(1, int(np.ceil((lat_max - lat_min) / target_resolution)))
            out_width = max(1, int(np.ceil((lon_max - lon_min) / target_resolution)))
        print(f'    Output grid: {out_width} x {out_height} (res={target_resolution:.6f} deg)')
        src_points = np.column_stack((lat_valid, lon_valid))
        tree = cKDTree(src_points)
        out_lats = np.linspace(lat_max - target_resolution / 2, lat_min + target_resolution / 2, out_height)
        out_lons = np.linspace(lon_min + target_resolution / 2, lon_max - target_resolution / 2, out_width)
        grid_lon, grid_lat = np.meshgrid(out_lons, out_lats)
        target_points = np.column_stack((grid_lat.ravel(), grid_lon.ravel()))
        dist, indices = tree.query(target_points)
        max_dist = target_resolution * 2.0
        out_bands = []
        for band_data in data_bands:
            flat_valid = band_data[valid]
            out_flat = np.full(out_height * out_width, np.nan, dtype=np.float32)
            mask = dist <= max_dist
            out_flat[mask] = flat_valid[indices[mask]]
            out_bands.append(out_flat.reshape(out_height, out_width))
        transform = rasterio.transform.from_bounds(lon_min, lat_min, lon_max, lat_max, out_width, out_height)
        return (out_bands, transform, out_height, out_width)
def read_scaled_band(ds, varname):
    """Read a variable and apply scale_factor + add_offset if present.
    Detects whether netCDF4 has already auto-applied scaling by checking
    the raw data dtype — integer means raw, float means pre-scaled."""
    raw = nc_read_var(ds, varname)
    scale = nc_get_attr(ds, varname, 'scale_factor', 1.0)
    offset = nc_get_attr(ds, varname, 'add_offset', 0.0)
    fill = nc_get_attr(ds, varname, '_FillValue', None)
    needs_scaling = np.issubdtype(raw.dtype, np.integer)
    if fill is not None:
        if needs_scaling:
            fill_mask = raw == fill
        else:
            fill_mask = np.isnan(raw) if np.issubdtype(raw.dtype, np.floating) else raw == fill
    else:
        fill_mask = np.zeros_like(raw, dtype=bool)
    data = np.float32(raw)
    if needs_scaling and (scale != 1.0 or offset != 0.0):
            data = data * float(scale) + float(offset)
    data[fill_mask] = np.nan
    vmin = nc_get_attr(ds, varname, 'valid_min', None)
    vmax = nc_get_attr(ds, varname, 'valid_max', None)
    if vmin is not None:
        ref = raw if needs_scaling else data
        data[ref < vmin] = np.nan
    if vmax is not None:
        ref = raw if needs_scaling else data
        data[ref > vmax] = np.nan
    return data
def read_geolocation(filepath, lat_var='latitude', lon_var='longitude'):
    """Read lat/lon arrays from a geolocation netCDF file."""
    ds = nc_open(filepath)
    lat = read_scaled_band(ds, lat_var)
    lon = read_scaled_band(ds, lon_var)
    nc_close(ds)
    return (lat, lon)
def build_envi_header(hdr_path, height, width, num_bands, dtype, crs, transform, wavelengths=None, band_names=None, description='', interleave='bsq'):
    """Write a custom ENVI header file."""
    dtype_map = {'uint8': 1, 'int16': 2, 'int32': 3, 'float32': 4, 'float64': 5, 'uint16': 12, 'uint32': 13}
    envi_dtype = dtype_map.get(np.dtype(dtype).name, 4)
    with open(hdr_path, 'w') as f:
        f.write('ENVI\n')
        f.write(f'description = {{{description}}}\n')
        f.write(f'samples = {width}\n')
        f.write(f'lines = {height}\n')
        f.write(f'bands = {num_bands}\n')
        f.write('header offset = 0\n')
        f.write('file type = ENVI Standard\n')
        f.write(f'data type = {envi_dtype}\n')
        f.write(f'interleave = {interleave}\n')
        f.write('byte order = 0\n')
        if transform is not None:
            x_size = abs(transform.a)
            y_size = abs(transform.e)
            ul_x = transform.c
            ul_y = transform.f
            if crs and crs.is_geographic:
                    f.write(f'map info = {{Geographic Lat/Lon, 1, 1, {ul_x}, {ul_y}, {x_size}, {y_size}, WGS-84}}\n')
        if crs is not None:
            f.write(f'coordinate system string = {{{crs.to_wkt()}}}\n')
        if wavelengths is not None:
            f.write('wavelength units = Micrometers\n')
            wl_str = ', '.join((f'{w:.4f}' for w in wavelengths))
            f.write(f'wavelength = {{\n  {wl_str}}}\n')
        if band_names is not None:
            bn_str = ',\n  '.join(band_names)
            f.write(f'band names = {{\n  {bn_str}}}\n')
        f.write('data ignore value = nan\n')
def write_multiband(band_data_list, output_path, transform, height, width, band_names, wavelengths=None, description='', write_geotiff=True, write_envi=True):
    """Write stacked bands as multiband GeoTIFF and/or ENVI in subfolders."""
    output_files = []
    num_bands = len(band_data_list)
    crs = build_wgs84_crs()
    out_dir = os.path.dirname(output_path)
    out_name = os.path.basename(output_path)
    stack = np.zeros((num_bands, height, width), dtype=np.float32)
    for i, data in enumerate(band_data_list):
        stack[i] = data
    if write_geotiff:
        gtiff_dir = os.path.join(out_dir, 'GeoTIFF')
        os.makedirs(gtiff_dir, exist_ok=True)
        gtiff_path = os.path.join(gtiff_dir, out_name + '.tif')
        profile = {'driver': 'GTiff', 'dtype': 'float32', 'width': width, 'height': height, 'count': num_bands, 'crs': crs, 'transform': transform, 'nodata': np.nan, 'compress': 'deflate', 'predictor': 2, 'zlevel': 6, 'tiled': True, 'blockxsize': 256, 'blockysize': 256}
        print(f'    Writing GeoTIFF: {os.path.basename(gtiff_path)}')
        with rasterio.open(gtiff_path, 'w', **profile) as dst:
            for i in range(num_bands):
                dst.write(stack[i], i + 1)
                dst.set_band_description(i + 1, band_names[i])
        output_files.append(gtiff_path)
    if write_envi:
        envi_dir = os.path.join(out_dir, 'ENVI')
        os.makedirs(envi_dir, exist_ok=True)
        envi_dat = os.path.join(envi_dir, out_name + '.dat')
        envi_hdr = os.path.join(envi_dir, out_name + '.hdr')
        print(f'    Writing ENVI: {os.path.basename(envi_dat)}')
        stack.tofile(envi_dat)
        build_envi_header(envi_hdr, height=height, width=width, num_bands=num_bands, dtype=np.float32, crs=crs, transform=transform, wavelengths=wavelengths, band_names=band_names, description=description, interleave='bsq')
        output_files.append(envi_dat)
        output_files.append(envi_hdr)
    return output_files
def extract_scene_id(folder_name):
    """
    Extract a compact scene ID from a Sentinel-3 .SEN3 folder name.
    e.g. S3A_SL_1_RBT____20260515T025525_... -> S3A_SL1RBT_20260515T025525
    """
    m = re.match('(S3[A-D])_([A-Z]{2})_(\\d)_([A-Z]{3})____(\\d{8}T\\d{6})', folder_name)
    if m:
        platform = m.group(1)
        inst = m.group(2)
        level = m.group(3)
        prod = m.group(4)
        timestamp = m.group(5)
        return f'{platform}_{inst}{level}{prod}_{timestamp}'
    else:
        return folder_name.replace('.SEN3', '')[:60]
def process_slstr_l1(scene_dir, output_dir, write_geotiff=True, write_envi=True, overwrite=True):
    """Process SLSTR L1B RBT: VNIR/SWIR radiance + TIR BT."""
    folder_name = os.path.basename(scene_dir.rstrip('/\\'))
    scene_id = extract_scene_id(folder_name)
    all_outputs = []
    config = SLSTR_VNIR_SWIR
    geo_path = os.path.join(scene_dir, config['geo_file'])
    if not os.path.isfile(geo_path):
        print(f"    WARNING: {config['geo_file']} not found, skipping VNIR/SWIR.")
    else:
        out_base = os.path.join(output_dir, f'{scene_id}_VNIR_SWIR_Radiance')
        if not overwrite and os.path.isfile(os.path.join(os.path.dirname(out_base), 'GeoTIFF', os.path.basename(out_base) + '.tif')):
            print('    VNIR/SWIR output exists, skipping.')
        else:
            print('\n  === SLSTR VNIR/SWIR Radiance (500m nadir) ===')
            lat, lon = read_geolocation(geo_path, 'latitude_an', 'longitude_an')
            print(f'    Geolocation: {lat.shape}')
            band_data = []
            actual_names = []
            actual_wl = []
            for i, band in enumerate(config['bands']):
                nc_file = os.path.join(scene_dir, f"{band}_radiance_{config['grid_suffix']}.nc")
                if not os.path.isfile(nc_file):
                    print(f'    WARNING: {band} not found, skipping.')
                    continue
                else:
                    varname = f"{band}_radiance_{config['grid_suffix']}"
                    print(f'    Reading {varname}...')
                    ds = nc_open(nc_file)
                    if nc_has_var(ds, varname):
                        data = read_scaled_band(ds, varname)
                        band_data.append(data)
                        actual_names.append(config['names'][i])
                        actual_wl.append(config['wavelengths'][i])
                    else:
                        print(f'    WARNING: Variable {varname} not found.')
                    nc_close(ds)
            if band_data:
                print(f'    Reprojecting {len(band_data)} bands...')
                out_bands, transform, oh, ow = reproject_swath_to_geographic(lat, lon, band_data)
                if out_bands:
                    outputs = write_multiband(out_bands, out_base, transform, oh, ow, band_names=actual_names, wavelengths=actual_wl, description=f'SLSTR L1B VNIR/SWIR Radiance ({scene_id})', write_geotiff=write_geotiff, write_envi=write_envi)
                    all_outputs.extend(outputs)
    config = SLSTR_TIR
    geo_path = os.path.join(scene_dir, config['geo_file'])
    if not os.path.isfile(geo_path):
        print(f"    WARNING: {config['geo_file']} not found, skipping TIR.")
    else:
        out_base = os.path.join(output_dir, f'{scene_id}_TIR_BT')
        if not overwrite and os.path.isfile(os.path.join(os.path.dirname(out_base), 'GeoTIFF', os.path.basename(out_base) + '.tif')):
            print('    TIR output exists, skipping.')
        else:
            print('\n  === SLSTR TIR Brightness Temperature (1km nadir) ===')
            lat, lon = read_geolocation(geo_path, 'latitude_in', 'longitude_in')
            print(f'    Geolocation: {lat.shape}')
            band_data = []
            actual_names = []
            actual_wl = []
            for i, band in enumerate(config['bands']):
                nc_file = os.path.join(scene_dir, f"{band}_BT_{config['grid_suffix']}.nc")
                if not os.path.isfile(nc_file):
                    print(f'    WARNING: {band}_BT not found, skipping.')
                    continue
                else:
                    varname = f"{band}_BT_{config['grid_suffix']}"
                    print(f'    Reading {varname}...')
                    ds = nc_open(nc_file)
                    if nc_has_var(ds, varname):
                        data = read_scaled_band(ds, varname)
                        band_data.append(data)
                        actual_names.append(config['names'][i])
                        actual_wl.append(config['wavelengths'][i])
                    else:
                        print(f'    WARNING: Variable {varname} not found.')
                    nc_close(ds)
            if band_data:
                print(f'    Reprojecting {len(band_data)} bands...')
                out_bands, transform, oh, ow = reproject_swath_to_geographic(lat, lon, band_data)
                if out_bands:
                    outputs = write_multiband(out_bands, out_base, transform, oh, ow, band_names=actual_names, wavelengths=actual_wl, description=f'SLSTR L1B TIR Brightness Temp ({scene_id})', write_geotiff=write_geotiff, write_envi=write_envi)
                    all_outputs.extend(outputs)
    config = SLSTR_FIRE
    geo_path = os.path.join(scene_dir, config['geo_file'])
    if not os.path.isfile(geo_path):
        print(f"    WARNING: {config['geo_file']} not found, skipping fire bands.")
    else:
        out_base = os.path.join(output_dir, f'{scene_id}_Fire_BT')
        if not overwrite and os.path.isfile(os.path.join(os.path.dirname(out_base), 'GeoTIFF', os.path.basename(out_base) + '.tif')):
            print('    Fire output exists, skipping.')
        else:
            print('\n  === SLSTR Fire BT (1km nadir, F1/F2) ===')
            lat, lon = read_geolocation(geo_path, 'latitude_fn', 'longitude_fn')
            print(f'    Geolocation: {lat.shape}')
            band_data = []
            actual_names = []
            actual_wl = []
            for i, band in enumerate(config['bands']):
                nc_file = os.path.join(scene_dir, f"{band}_BT_{config['grid_suffix']}.nc")
                if not os.path.isfile(nc_file):
                    print(f'    WARNING: {band}_BT not found, skipping.')
                    continue
                else:
                    varname = f"{band}_BT_{config['grid_suffix']}"
                    print(f'    Reading {varname}...')
                    ds = nc_open(nc_file)
                    if nc_has_var(ds, varname):
                        data = read_scaled_band(ds, varname)
                        band_data.append(data)
                        actual_names.append(config['names'][i])
                        actual_wl.append(config['wavelengths'][i])
                    else:
                        print(f'    WARNING: Variable {varname} not found.')
                    nc_close(ds)
            if band_data:
                print(f'    Reprojecting {len(band_data)} bands...')
                out_bands, transform, oh, ow = reproject_swath_to_geographic(lat, lon, band_data)
                if out_bands:
                    outputs = write_multiband(out_bands, out_base, transform, oh, ow, band_names=actual_names, wavelengths=actual_wl, description=f'SLSTR L1B Fire BT ({scene_id})', write_geotiff=write_geotiff, write_envi=write_envi)
                    all_outputs.extend(outputs)
    return all_outputs
def process_slstr_l2_lst(scene_dir, output_dir, write_geotiff=True, write_envi=True, overwrite=True):
    """Process SLSTR L2 LST."""
    folder_name = os.path.basename(scene_dir.rstrip('/\\'))
    scene_id = extract_scene_id(folder_name)
    all_outputs = []
    out_base = os.path.join(output_dir, f'{scene_id}_LST')
    if not overwrite and os.path.isfile(os.path.join(os.path.dirname(out_base), 'GeoTIFF', os.path.basename(out_base) + '.tif')):
        print('    Output exists, skipping.')
        return all_outputs
    else:
        print('\n  === SLSTR L2 Land Surface Temperature (1km) ===')
        geo_path = os.path.join(scene_dir, 'geodetic_in.nc')
        if not os.path.isfile(geo_path):
            print('    ERROR: geodetic_in.nc not found.')
            return all_outputs
        else:
            lat, lon = read_geolocation(geo_path, 'latitude_in', 'longitude_in')
            print(f'    Geolocation: {lat.shape}')
            lst_path = os.path.join(scene_dir, 'LST_in.nc')
            if not os.path.isfile(lst_path):
                print('    ERROR: LST_in.nc not found.')
                return all_outputs
            else:
                band_data = []
                band_names = []
                band_wl = []
                ds = nc_open(lst_path)
                if nc_has_var(ds, 'LST'):
                    print('    Reading LST...')
                    lst = read_scaled_band(ds, 'LST')
                    band_data.append(lst)
                    band_names.append('Land Surface Temperature (K)')
                    band_wl.append(11.0)
                if nc_has_var(ds, 'LST_uncertainty'):
                    print('    Reading LST_uncertainty...')
                    unc = read_scaled_band(ds, 'LST_uncertainty')
                    band_data.append(unc)
                    band_names.append('LST Uncertainty (K)')
                    band_wl.append(11.0)
                nc_close(ds)
                if band_data:
                    print(f'    Reprojecting {len(band_data)} bands...')
                    out_bands, transform, oh, ow = reproject_swath_to_geographic(lat, lon, band_data)
                    if out_bands:
                        outputs = write_multiband(out_bands, out_base, transform, oh, ow, band_names=band_names, wavelengths=band_wl, description=f'SLSTR L2 LST ({scene_id})', write_geotiff=write_geotiff, write_envi=write_envi)
                        all_outputs.extend(outputs)
                return all_outputs
def process_olci_l1(scene_dir, output_dir, write_geotiff=True, write_envi=True, overwrite=True):
    """Process OLCI L1B EFR: 21-band TOA radiance."""
    folder_name = os.path.basename(scene_dir.rstrip('/\\'))
    scene_id = extract_scene_id(folder_name)
    all_outputs = []
    out_base = os.path.join(output_dir, f'{scene_id}_TOA_Radiance')
    if not overwrite and os.path.isfile(os.path.join(os.path.dirname(out_base), 'GeoTIFF', os.path.basename(out_base) + '.tif')):
        print('    Output exists, skipping.')
        return all_outputs
    else:
        print('\n  === OLCI L1B TOA Radiance (300m, 21 bands) ===')
        config = OLCI_BANDS
        geo_path = os.path.join(scene_dir, config['geo_file'])
        if not os.path.isfile(geo_path):
            print(f"    ERROR: {config['geo_file']} not found.")
            return all_outputs
        else:
            lat, lon = read_geolocation(geo_path)
            print(f'    Geolocation: {lat.shape}')
            band_data = []
            actual_names = []
            actual_wl = []
            for i, band in enumerate(config['bands']):
                nc_file = os.path.join(scene_dir, f'{band}_radiance.nc')
                if not os.path.isfile(nc_file):
                    print(f'    WARNING: {band} not found, skipping.')
                    continue
                else:
                    varname = f'{band}_radiance'
                    print(f'    Reading {varname}...')
                    ds = nc_open(nc_file)
                    if nc_has_var(ds, varname):
                        data = read_scaled_band(ds, varname)
                        band_data.append(data)
                        actual_names.append(config['names'][i])
                        actual_wl.append(config['wavelengths'][i])
                    else:
                        print(f'    WARNING: Variable {varname} not found.')
                    nc_close(ds)
            if band_data:
                print(f'    Reprojecting {len(band_data)} bands...')
                out_bands, transform, oh, ow = reproject_swath_to_geographic(lat, lon, band_data)
                if out_bands:
                    outputs = write_multiband(out_bands, out_base, transform, oh, ow, band_names=actual_names, wavelengths=actual_wl, description=f'OLCI L1B TOA Radiance ({scene_id})', write_geotiff=write_geotiff, write_envi=write_envi)
                    all_outputs.extend(outputs)
            return all_outputs
def process_synergy_l2(scene_dir, output_dir, write_geotiff=True, write_envi=True, overwrite=True):
    """Process Synergy L2 SYN: OLCI + SLSTR surface reflectance."""
    folder_name = os.path.basename(scene_dir.rstrip('/\\'))
    scene_id = extract_scene_id(folder_name)
    all_outputs = []
    config = SYN_OLCI_BANDS
    out_base = os.path.join(output_dir, f'{scene_id}_OLCI_SurfRefl')
    if not overwrite and os.path.isfile(os.path.join(os.path.dirname(out_base), 'GeoTIFF', os.path.basename(out_base) + '.tif')):
        print('    OLCI reflectance output exists, skipping.')
    else:
        print('\n  === Synergy L2 OLCI Surface Reflectance (300m) ===')
        geo_path = os.path.join(scene_dir, config['geo_file'])
        if not os.path.isfile(geo_path):
            print(f"    ERROR: {config['geo_file']} not found.")
        else:
            lat, lon = read_geolocation(geo_path, 'lat', 'lon')
            print(f'    Geolocation: {lat.shape}')
            band_data = []
            actual_names = []
            actual_wl = []
            for i, band in enumerate(config['bands']):
                nc_file = os.path.join(scene_dir, f'{band}_reflectance.nc')
                if not os.path.isfile(nc_file):
                    print(f'    WARNING: {band} not found, skipping.')
                    continue
                else:
                    varname = band.replace('Syn_', 'SDR_')
                    print(f'    Reading {varname}...')
                    ds = nc_open(nc_file)
                    if nc_has_var(ds, varname):
                        data = read_scaled_band(ds, varname)
                        band_data.append(data)
                        actual_names.append(config['names'][i])
                        actual_wl.append(config['wavelengths'][i])
                    else:
                        print(f'    WARNING: Variable {varname} not found.')
                    nc_close(ds)
            if band_data:
                print(f'    Reprojecting {len(band_data)} bands...')
                out_bands, transform, oh, ow = reproject_swath_to_geographic(lat, lon, band_data)
                if out_bands:
                    outputs = write_multiband(out_bands, out_base, transform, oh, ow, band_names=actual_names, wavelengths=actual_wl, description=f'Synergy L2 OLCI Surface Reflectance ({scene_id})', write_geotiff=write_geotiff, write_envi=write_envi)
                    all_outputs.extend(outputs)
    config = SYN_SLSTR_BANDS
    out_base = os.path.join(output_dir, f'{scene_id}_SLSTR_SurfRefl')
    if not overwrite and os.path.isfile(os.path.join(os.path.dirname(out_base), 'GeoTIFF', os.path.basename(out_base) + '.tif')):
        print('    SLSTR reflectance output exists, skipping.')
    else:
        print('\n  === Synergy L2 SLSTR Surface Reflectance (300m) ===')
        geo_path = os.path.join(scene_dir, config['geo_file'])
        if os.path.isfile(geo_path):
            lat, lon = read_geolocation(geo_path, 'lat', 'lon')
            band_data = []
            actual_names = []
            actual_wl = []
            for i, band in enumerate(config['bands']):
                nc_file = os.path.join(scene_dir, f'{band}_reflectance.nc')
                if not os.path.isfile(nc_file):
                    print(f'    WARNING: {band} not found, skipping.')
                    continue
                else:
                    varname = band.replace('Syn_', 'SDR_')
                    print(f'    Reading {varname}...')
                    ds = nc_open(nc_file)
                    if nc_has_var(ds, varname):
                        data = read_scaled_band(ds, varname)
                        band_data.append(data)
                        actual_names.append(config['names'][i])
                        actual_wl.append(config['wavelengths'][i])
                    else:
                        print(f'    WARNING: Variable {varname} not found.')
                    nc_close(ds)
            if band_data:
                print(f'    Reprojecting {len(band_data)} bands...')
                out_bands, transform, oh, ow = reproject_swath_to_geographic(lat, lon, band_data)
                if out_bands:
                    outputs = write_multiband(out_bands, out_base, transform, oh, ow, band_names=actual_names, wavelengths=actual_wl, description=f'Synergy L2 SLSTR Surface Reflectance ({scene_id})', write_geotiff=write_geotiff, write_envi=write_envi)
                    all_outputs.extend(outputs)
    return all_outputs
def detect_product_type(folder_name):
    """Detect product type from .SEN3 folder name."""
    for product_type, pattern in PRODUCT_PATTERNS.items():
        if pattern.search(folder_name):
            return product_type
def discover_scenes(input_dir, product_filter=None):
    """
    Scan for Sentinel-3 .SEN3 folders.
    Also checks if input_dir itself is a .SEN3 folder.
    If product_filter is set (e.g. \'SL_1_RBT\'), only matching scenes
    are returned.
    Returns list of (scene_dir, product_type, folder_name).
    """
    scenes = []
    folder_name = os.path.basename(input_dir.rstrip('/\\'))
    product_type = detect_product_type(folder_name)
    if product_type:
        if product_filter is None or product_type == product_filter:
            scenes.append((input_dir, product_type, folder_name))
        return scenes
    else:
        if os.path.isfile(os.path.join(input_dir, 'xfdumanifest.xml')):
            parent_name = os.path.basename(input_dir.rstrip('/\\'))
            product_type = detect_product_type(parent_name)
            if product_type:
                if product_filter is None or product_type == product_filter:
                    scenes.append((input_dir, product_type, parent_name))
                return scenes
        for entry in sorted(os.listdir(input_dir)):
            subdir = os.path.join(input_dir, entry)
            if os.path.isdir(subdir):
                product_type = detect_product_type(entry)
                if product_type and (product_filter is None or product_type == product_filter):
                        scenes.append((subdir, product_type, entry))
        return scenes
def process_scene(scene_dir, product_type, output_dir, write_geotiff=True, write_envi=True, overwrite=True):
    """Process a single Sentinel-3 scene by product type."""
    if product_type == 'SL_1_RBT':
        return process_slstr_l1(scene_dir, output_dir, write_geotiff, write_envi, overwrite)
    else:
        if product_type == 'SL_2_LST':
            return process_slstr_l2_lst(scene_dir, output_dir, write_geotiff, write_envi, overwrite)
        else:
            if product_type == 'OL_1_EFR':
                return process_olci_l1(scene_dir, output_dir, write_geotiff, write_envi, overwrite)
            else:
                if product_type == 'SY_2_SYN':
                    return process_synergy_l2(scene_dir, output_dir, write_geotiff, write_envi, overwrite)
                else:
                    print(f'  WARNING: Unknown product type \'{product_type}\'.')
                    return []
def process_directory(input_dir, output_dir=None, write_geotiff=True, write_envi=True, overwrite=True, product_filter=None):
    """Process Sentinel-3 scenes. If product_filter is set, only matching
    product types are processed (e.g. \'SL_1_RBT\', \'OL_1_EFR\')."""
    if output_dir is None:
        output_dir = os.path.join(input_dir, 'converted')
    filter_label = f' (filter: {product_filter})' if product_filter else ''
    print(f'Scanning for Sentinel-3 scenes in: {input_dir}{filter_label}')
    scenes = discover_scenes(input_dir, product_filter=product_filter)
    if not scenes:
        print('No Sentinel-3 .SEN3 scenes found.')
        print('Supported: SL_1_RBT, SL_2_LST, OL_1_EFR, SY_2_SYN')
        return
    else:
        print(f'Found {len(scenes)} scene(s):\n')
        os.makedirs(output_dir, exist_ok=True)
        all_outputs = []
        for i, (scene_dir, product_type, folder_name) in enumerate(scenes, 1):
            print(f"{'======================================================================'}")
            print(f'Scene {i}/{len(scenes)}: {folder_name}')
            print(f'  Product: {product_type}')
            print(f"{'======================================================================'}")
            try:
                outputs = process_scene(scene_dir, product_type, output_dir, write_geotiff=write_geotiff, write_envi=write_envi, overwrite=overwrite)
                all_outputs.extend(outputs)
            except Exception as e:
                print(f'  ERROR: {e}')
                import traceback
                traceback.print_exc()
            print()
        print(f'Processing complete. {len(all_outputs)} output file(s) in: {output_dir}')
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Sentinel-3 Product Converter')
        print('==================================================')
        print(f'\nUsage: python {os.path.basename(__file__)} <input_dir> [output_dir] [options]')
        print('\nSupported products (auto-detected from .SEN3 folder name):')
        print('  SL_1_RBT   SLSTR L1B Radiance + BT (S1-S9, F1-F2)')
        print('  SL_2_LST   SLSTR L2 Land Surface Temperature')
        print('  OL_1_EFR   OLCI L1B Full Res TOA Radiance (21 bands)')
        print('  SY_2_SYN   Synergy L2 Surface Reflectance (OLCI+SLSTR)')
        print('\nPlatforms: S3A, S3B, S3C, S3D')
        print('\nOptions: --geotiff, --envi (default: both)')
        sys.exit(0)
    args = sys.argv[1:]
    write_geotiff = True
    write_envi = True
    positional = []
    for arg in args:
        if arg.lower() == '--geotiff':
            write_geotiff, write_envi = (True, False)
        else:
            if arg.lower() == '--envi':
                write_geotiff, write_envi = (False, True)
            else:
                positional.append(arg)
    input_dir = positional[0]
    output_dir = positional[1] if len(positional) > 1 else None
    if not os.path.isdir(input_dir):
        print(f'Error: {input_dir} does not exist')
        sys.exit(1)
    process_directory(input_dir, output_dir, write_geotiff=write_geotiff, write_envi=write_envi)