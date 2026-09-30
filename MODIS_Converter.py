# Decompiled with PyLingual (https://pylingual.io)
# Internal filename: 'MODIS_Converter.py'
# Bytecode version: 3.11a7e (3495)
# Source timestamp: 1970-01-01 00:00:00 UTC (0)

"""
MODIS Product Converter
=========================
Reads MODIS HDF-EOS2 gridded products and outputs multiband GeoTIFF and/or
ENVI files in the native sinusoidal projection.

Supported products:
  MOD09GA / MYD09GA (Surface Reflectance, daily, 500m)
    -> 7-band reflectance (Bands 1-7) at 500m
  MOD11A1 / MYD11A1 (LST/Emissivity, daily, 1km)
    -> 4-band file: LST_Day, LST_Night, Emis_31, Emis_32 at 1km

Auto-detects product type from filename.
Supports both Terra (MOD) and Aqua (MYD) variants.

Output is kept in native sinusoidal projection — ArcGIS Pro / QGIS
can reproject as needed.

Dependencies: numpy, rasterio, pyhdf
"""
import os
import re
import sys
import numpy as np
try:
    import rasterio
    from rasterio.crs import CRS
    from rasterio.transform import from_bounds
except ImportError:
    print('ERROR: rasterio is required. Install with: pip install rasterio')
    sys.exit(1)
try:
    from pyhdf.SD import SD, SDC
except ImportError:
    print('ERROR: pyhdf is required. Install with: pip install pyhdf')
    sys.exit(1)
MODIS_SINUSOIDAL_WKT = 'PROJCS[\"MODIS Sinusoidal\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6371007.181,0]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Sinusoidal\"],PARAMETER[\"central_meridian\",0],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"Meter\",1]]'
def build_sinusoidal_crs():
    """Build the MODIS sinusoidal CRS from WKT (avoids PROJ database)."""
    return CRS.from_wkt(MODIS_SINUSOIDAL_WKT)
MODIS_FILE_PATTERN = re.compile('^(MOD|MYD)(09GA|11A1)\\.A(\\d{7})\\.h(\\d{2})v(\\d{2})\\.\\d{3}\\.\\d+\\.hdf$', re.IGNORECASE)
MOD09_BANDS = {'datasets': ['sur_refl_b01_1', 'sur_refl_b02_1', 'sur_refl_b03_1', 'sur_refl_b04_1', 'sur_refl_b05_1', 'sur_refl_b06_1', 'sur_refl_b07_1'], 'scale_factor': 0.0001, 'add_offset': 0.0, 'fill_value': (-28672), 'valid_range': ((-100), 16000), 'grid_name': 'MODIS_Grid_500m_2D', 'wavelengths': [0.645, 0.858, 0.469, 0.555, 1.24, 1.64, 2.13], 'band_names': ['Band 1 - Red (0.645 um)', 'Band 2 - NIR (0.858 um)', 'Band 3 - Blue (0.469 um)', 'Band 4 - Green (0.555 um)', 'Band 5 - SWIR-1 (1.240 um)', 'Band 6 - SWIR-2 (1.640 um)', 'Band 7 - SWIR-3 (2.130 um)'], 'output_suffix': 'SurfRefl', 'description': 'MODIS Surface Reflectance (Bands 1-7)'}
MOD11_BANDS = {'datasets': ['LST_Day_1km', 'LST_Night_1km', 'Emis_31', 'Emis_32'], 'scale_factors': [0.02, 0.02, 0.002, 0.002], 'add_offsets': [0.0, 0.0, 0.49, 0.49], 'fill_values': [0, 0, 0, 0], 'valid_ranges': [(7500, 65535), (7500, 65535), (1, 255), (1, 255)], 'grid_name': 'MODIS_Grid_Daily_1km_LST', 'wavelengths': [11.03, 11.03, 11.03, 12.02], 'band_names': ['LST Day (K)', 'LST Night (K)', 'Emissivity Band 31 (11.03 um)', 'Emissivity Band 32 (12.02 um)'], 'output_suffix': 'LST_Emis', 'description': 'MODIS Land Surface Temperature and Emissivity'}
def parse_struct_metadata(hdf_file):
    """
    Parse StructMetadata.0 to extract grid definitions.
    Returns dict of grid_name -> {ul_x, ul_y, lr_x, lr_y, xdim, ydim}
    """
    sd = SD(hdf_file, SDC.READ)
    attrs = sd.attributes()
    struct_meta = attrs.get('StructMetadata.0', '')
    sd.end()
    grids = {}
    grid_pattern = re.compile('GROUP=GRID_\\d+.*?GridName=\"([^\"]+)\".*?XDim=(\\d+).*?YDim=(\\d+).*?UpperLeftPointMtrs=\\(([-\\d.]+),([-\\d.]+)\\).*?LowerRightMtrs=\\(([-\\d.]+),([-\\d.]+)\\)', re.DOTALL)
    for match in grid_pattern.finditer(struct_meta):
        name = match.group(1)
        grids[name] = {'xdim': int(match.group(2)), 'ydim': int(match.group(3)), 'ul_x': float(match.group(4)), 'ul_y': float(match.group(5)), 'lr_x': float(match.group(6)), 'lr_y': float(match.group(7))}
    return grids
def build_transform(grid_info):
    """Build a rasterio affine transform from grid corner coordinates."""
    ul_x = grid_info['ul_x']
    ul_y = grid_info['ul_y']
    lr_x = grid_info['lr_x']
    lr_y = grid_info['lr_y']
    xdim = grid_info['xdim']
    ydim = grid_info['ydim']
    pixel_x = (lr_x - ul_x) / xdim
    pixel_y = (lr_y - ul_y) / ydim
    return rasterio.transform.Affine(pixel_x, 0, ul_x, 0, pixel_y, ul_y)
def read_hdf_dataset(hdf_path, dataset_name):
    """Read a single dataset from an HDF-EOS2 file."""
    sd = SD(hdf_path, SDC.READ)
    try:
        ds = sd.select(dataset_name)
        data = ds[:].astype(np.float64)
        attrs = ds.attributes()
        ds.endaccess()
    except Exception:
        sd.end()
        return (None, {})
    sd.end()
    return (data, attrs)
def apply_scaling(data, scale_factor, add_offset, fill_value, valid_range):
    """
    Apply scale/offset and mask invalid pixels.
    physical = DN * scale_factor + add_offset
    """
    result = np.float32(data)
    fill_mask = data == fill_value
    if valid_range is not None:
        vmin, vmax = valid_range
        range_mask = (data < vmin) | (data > vmax)
        fill_mask |= range_mask
    result = result * scale_factor + add_offset
    result[fill_mask] = np.nan
    return result
def discover_files(input_dir):
    """
    Scan for MODIS HDF files. Returns list of (filepath, product_type, info)
    where product_type is \'MOD09GA\', \'MYD09GA\', \'MOD11A1\', or \'MYD11A1\'.
    """
    files = []
    for fname in sorted(os.listdir(input_dir)):
        match = MODIS_FILE_PATTERN.match(fname)
        if match:
            platform = match.group(1).upper()
            product = match.group(2).upper()
            product_type = f'{platform}{product}'
            filepath = os.path.join(input_dir, fname)
            files.append((filepath, product_type, fname))
    return files
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
            f.write(f'map info = {{Sinusoidal, 1, 1, {ul_x}, {ul_y}, {x_size}, {y_size}, units=Meters}}\n')
            # ENVI needs explicit projection parameters for non-UTM projections.
            # Format: {type, sphere radius, central meridian, false E, false N, name, units}
            # Type 16 = Sinusoidal; radius is the MODIS authalic sphere.
            f.write('projection info = {16, 6371007.181, 0.0, 0.0, 0.0, '
                    'Sinusoidal, units=Meters}\n')
        # Deliberately no 'coordinate system string' here. GDAL normalizes the
        # MODIS sinusoidal WKT to SPHEROID["WGS 84", 6371007.181, 0,
        # AUTHORITY["EPSG","7030"]] — the inline radius is the correct authalic
        # sphere, but EPSG:7030 is the real WGS 84 ellipsoid (6378137,
        # 1/298.257). ENVI resolves the AUTHORITY code in preference to the
        # inline value and so inverts the projection on an ellipsoid instead of
        # a sphere, putting the lat/lon readout wildly out. ArcGIS tolerates the
        # contradiction, which is why GeoTIFF output was unaffected.
        # map info + projection info above fully define the grid for ENVI.
        if wavelengths is not None:
            f.write('wavelength units = Micrometers\n')
            wl_str = ', '.join((f'{w:.3f}' for w in wavelengths))
            f.write(f'wavelength = {{\n  {wl_str}}}\n')
        if band_names is not None:
            bn_str = ',\n  '.join(band_names)
            f.write(f'band names = {{\n  {bn_str}}}\n')
        f.write('data ignore value = nan\n')
def write_multiband(band_data_list, output_path, height, width, crs, transform, band_names, wavelengths=None, description='', write_geotiff=True, write_envi=True):
    """Write stacked bands as multiband GeoTIFF and/or ENVI."""
    output_files = []
    num_bands = len(band_data_list)
    stack = np.zeros((num_bands, height, width), dtype=np.float32)
    for i, data in enumerate(band_data_list):
        stack[i] = data
    if write_geotiff:
        gtiff_path = output_path + '.tif'
        profile = {'driver': 'GTiff', 'dtype': 'float32', 'width': width, 'height': height, 'count': num_bands, 'crs': crs, 'transform': transform, 'nodata': np.nan, 'compress': 'deflate', 'predictor': 2, 'zlevel': 6, 'tiled': True, 'blockxsize': 256, 'blockysize': 256}
        print(f'    Writing GeoTIFF: {os.path.basename(gtiff_path)}')
        with rasterio.open(gtiff_path, 'w', **profile) as dst:
            for i in range(num_bands):
                dst.write(stack[i], i + 1)
                dst.set_band_description(i + 1, band_names[i])
        output_files.append(gtiff_path)
    if write_envi:
        envi_dat = output_path + '.dat'
        envi_hdr = output_path + '.hdr'
        print(f'    Writing ENVI: {os.path.basename(envi_dat)}')
        stack.tofile(envi_dat)
        build_envi_header(envi_hdr, height=height, width=width, num_bands=num_bands, dtype=np.float32, crs=crs, transform=transform, wavelengths=wavelengths, band_names=band_names, description=description, interleave='bsq')
        output_files.append(envi_dat)
        output_files.append(envi_hdr)
    return output_files
def process_mod09(hdf_path, output_dir, write_geotiff=True, write_envi=True, overwrite=True):
    """Process a MOD09GA or MYD09GA file (Surface Reflectance)."""
    config = MOD09_BANDS
    fname = os.path.basename(hdf_path)
    scene_id = fname.replace('.hdf', '')
    out_base = os.path.join(output_dir, f"{scene_id}_{config['output_suffix']}")
    if not overwrite and os.path.isfile(out_base + '.tif'):
        print('    Output exists, skipping (overwrite=False)')
        return []
    else:
        grids = parse_struct_metadata(hdf_path)
        grid_info = grids.get(config['grid_name'])
        if grid_info is None:
            print(f"    ERROR: Grid \'{config['grid_name']}\' not found in StructMetadata.")
            return []
        else:
            transform = build_transform(grid_info)
            crs = build_sinusoidal_crs()
            height = grid_info['ydim']
            width = grid_info['xdim']
            print(f'    Grid: {width} x {height}, pixel: {abs(transform.a):.1f} x {abs(transform.e):.1f} m')
            band_data = []
            actual_names = []
            actual_wl = []
            for i, ds_name in enumerate(config['datasets']):
                data, attrs = read_hdf_dataset(hdf_path, ds_name)
                if data is None:
                    print(f'    WARNING: Dataset \'{ds_name}\' not found, skipping.')
                    continue
                else:
                    sf = attrs.get('scale_factor', config['scale_factor'])
                    ao = attrs.get('add_offset', config['add_offset'])
                    fv = attrs.get('_FillValue', config['fill_value'])
                    vr = attrs.get('valid_range', config['valid_range'])
                    if isinstance(vr, np.ndarray):
                        vr = (int(vr[0]), int(vr[1]))
                    print(f'    Reading {ds_name} (scale={sf}, offset={ao})...')
                    scaled = apply_scaling(data, sf, ao, fv, vr)
                    band_data.append(scaled)
                    actual_names.append(config['band_names'][i])
                    actual_wl.append(config['wavelengths'][i])
            if not band_data:
                print('    ERROR: No bands successfully read.')
                return []
            else:
                return write_multiband(band_data, out_base, height, width, crs, transform, band_names=actual_names, wavelengths=actual_wl, description=f"{config['description']} ({scene_id})", write_geotiff=write_geotiff, write_envi=write_envi)
def process_mod11(hdf_path, output_dir, write_geotiff=True, write_envi=True, overwrite=True):
    """Process a MOD11A1 or MYD11A1 file (LST + Emissivity)."""
    config = MOD11_BANDS
    fname = os.path.basename(hdf_path)
    scene_id = fname.replace('.hdf', '')
    out_base = os.path.join(output_dir, f"{scene_id}_{config['output_suffix']}")
    if not overwrite and os.path.isfile(out_base + '.tif'):
        print('    Output exists, skipping (overwrite=False)')
        return []
    else:
        grids = parse_struct_metadata(hdf_path)
        grid_info = grids.get(config['grid_name'])
        if grid_info is None:
            print(f"    ERROR: Grid \'{config['grid_name']}\' not found in StructMetadata.")
            return []
        else:
            transform = build_transform(grid_info)
            crs = build_sinusoidal_crs()
            height = grid_info['ydim']
            width = grid_info['xdim']
            print(f'    Grid: {width} x {height}, pixel: {abs(transform.a):.1f} x {abs(transform.e):.1f} m')
            band_data = []
            actual_names = []
            actual_wl = []
            for i, ds_name in enumerate(config['datasets']):
                data, attrs = read_hdf_dataset(hdf_path, ds_name)
                if data is None:
                    print(f'    WARNING: Dataset \'{ds_name}\' not found, skipping.')
                    continue
                else:
                    sf = config['scale_factors'][i]
                    ao = config['add_offsets'][i]
                    fv = config['fill_values'][i]
                    vr = config['valid_ranges'][i]
                    print(f'    Reading {ds_name} (scale={sf}, offset={ao})...')
                    scaled = apply_scaling(data, sf, ao, fv, vr)
                    band_data.append(scaled)
                    actual_names.append(config['band_names'][i])
                    actual_wl.append(config['wavelengths'][i])
            if not band_data:
                print('    ERROR: No bands successfully read.')
                return []
            else:
                return write_multiband(band_data, out_base, height, width, crs, transform, band_names=actual_names, wavelengths=actual_wl, description=f"{config['description']} ({scene_id})", write_geotiff=write_geotiff, write_envi=write_envi)
def process_file(hdf_path, output_dir, write_geotiff=True, write_envi=True, overwrite=True):
    """Process a single MODIS HDF file (auto-detects product type)."""
    fname = os.path.basename(hdf_path)
    match = MODIS_FILE_PATTERN.match(fname)
    if not match:
        print(f'  WARNING: \'{fname}\' does not match MODIS filename pattern.')
        return []
    else:
        product = match.group(2).upper()
        if product == '09GA':
            return process_mod09(hdf_path, output_dir, write_geotiff=write_geotiff, write_envi=write_envi, overwrite=overwrite)
        else:
            if product == '11A1':
                return process_mod11(hdf_path, output_dir, write_geotiff=write_geotiff, write_envi=write_envi, overwrite=overwrite)
            else:
                print(f'  WARNING: Unsupported MODIS product \'{product}\'.')
                return []
def process_directory(input_dir, output_dir=None, write_geotiff=True, write_envi=True, overwrite=True):
    """Process all MODIS HDF files in a directory."""
    if output_dir is None:
        output_dir = os.path.join(input_dir, 'converted')
    print(f'Scanning for MODIS files in: {input_dir}')
    files = discover_files(input_dir)
    if not files:
        print('No MODIS files found.')
        print('Supported: MOD09GA, MYD09GA, MOD11A1, MYD11A1 (.hdf)')
        return
    else:
        print(f'Found {len(files)} file(s):\n')
        os.makedirs(output_dir, exist_ok=True)
        all_outputs = []
        for i, (filepath, product_type, fname) in enumerate(files, 1):
            print(f"{'======================================================================'}")
            print(f'File {i}/{len(files)}: {fname}')
            print(f'  Product: {product_type}')
            print(f"{'======================================================================'}")
            try:
                outputs = process_file(filepath, output_dir, write_geotiff=write_geotiff, write_envi=write_envi, overwrite=overwrite)
                all_outputs.extend(outputs)
            except Exception as e:
                print(f'  ERROR processing {fname}: {e}')
                import traceback
                traceback.print_exc()
            print()
        print(f'Processing complete. {len(all_outputs)} output file(s) in: {output_dir}')
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('MODIS Product Converter')
        print('==================================================')
        print(f'\nUsage: python {os.path.basename(__file__)} <input_dir> [output_dir] [options]')
        print('\nSupported products (Terra + Aqua):')
        print('  MOD09GA / MYD09GA  Surface Reflectance (7 bands, 500m)')
        print('  MOD11A1 / MYD11A1  LST + Emissivity (4 bands, 1km)')
        print('\nOutput is in native sinusoidal projection.')
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