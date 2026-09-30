"""
ACOLITE L2W Converter
=====================
Converts ACOLITE L2W NetCDF files to multiband GeoTIFF and/or ENVI format.

Extracts surface reflectance (rhos_*) and/or TOA reflectance (rhot_*) bands,
stacks them in wavelength order, and writes georeferenced outputs.

Supports both Landsat and Sentinel-2 ACOLITE outputs.

Usage:
    python ACOLITE_L2W_Converter.py <input_file_or_folder> [output_folder] [options]

Options:
    --geotiff    GeoTIFF only (default: both)
    --envi       ENVI only
    --rhot       Extract TOA reflectance instead of surface reflectance
    --both       Extract both rhos and rhot as separate files
"""

import os
import re
import sys
import glob
import numpy as np

try:
    import h5py
    _HAS_H5PY = True
except ImportError:
    _HAS_H5PY = False

try:
    import rasterio
    from rasterio.crs import CRS
    from rasterio.transform import Affine
    _HAS_RASTERIO = True
except ImportError:
    _HAS_RASTERIO = False


def read_acolite_l2w(filepath, product='rhos'):
    """
    Read an ACOLITE L2W NetCDF file.

    Parameters
    ----------
    filepath : str
        Path to .nc file
    product : str
        'rhos' for surface reflectance, 'rhot' for TOA reflectance

    Returns
    -------
    dict with keys: arrays, band_names, wavelengths, transform, crs, height, width
    """
    prefix = f'{product}_'

    with h5py.File(filepath, 'r') as f:
        # Find all matching datasets
        datasets = sorted([k for k in f.keys() if k.startswith(prefix)
                           and k != f'{product}_ds'],
                          key=lambda k: float(re.findall(r'[\d.]+', k.split('_')[1])[0]))

        if not datasets:
            print(f"  No {prefix}* datasets found in {os.path.basename(filepath)}")
            return None

        # Read geolocation
        transform = None
        crs = None

        # Try to get projection info from attributes
        attrs = dict(f.attrs)

        # Look for projection key and scene range
        proj_key = attrs.get('proj4_string', None)
        if proj_key and isinstance(proj_key, bytes):
            proj_key = proj_key.decode('utf-8')

        # Try CRS from proj4 string
        if proj_key:
            try:
                crs = CRS.from_proj4(proj_key)
            except Exception:
                pass

        # Try CRS from crs_wkt attribute
        if crs is None:
            crs_wkt = attrs.get('crs_wkt', None)
            if crs_wkt:
                if isinstance(crs_wkt, bytes):
                    crs_wkt = crs_wkt.decode('utf-8')
                try:
                    crs = CRS.from_wkt(crs_wkt)
                except Exception:
                    pass

        # Try from transverse_mercator grid mapping
        if crs is None and 'transverse_mercator' in f:
            tm = f['transverse_mercator']
            tm_attrs = dict(tm.attrs)
            crs_wkt = tm_attrs.get('crs_wkt', None)
            if crs_wkt:
                if isinstance(crs_wkt, bytes):
                    crs_wkt = crs_wkt.decode('utf-8')
                try:
                    crs = CRS.from_wkt(crs_wkt)
                except Exception:
                    pass
            if crs is None:
                # Try building from zone
                zone = attrs.get('scene_zone', None)
                if zone:
                    zone = int(zone)
                    epsg = 32600 + zone
                    crs = CRS.from_epsg(epsg)

        # Build transform from x/y coordinates
        if 'x' in f and 'y' in f:
            x = f['x'][:]
            y = f['y'][:]
            if len(x) > 1 and len(y) > 1:
                dx = float(x[1] - x[0])
                dy = float(y[1] - y[0])  # typically negative
                transform = Affine(dx, 0, float(x[0]) - dx/2,
                                   0, dy, float(y[0]) - dy/2)
        elif 'lon' in f and 'lat' in f:
            lon = f['lon'][:]
            lat = f['lat'][:]
            if lon.ndim == 1 and len(lon) > 1:
                dx = float(lon[1] - lon[0])
                dy = float(lat[1] - lat[0])
                transform = Affine(dx, 0, float(lon[0]) - dx/2,
                                   0, dy, float(lat[0]) - dy/2)
                if crs is None:
                    crs = CRS.from_epsg(4326)

        # Read bands
        arrays = []
        band_names = []
        wavelengths = []

        for ds_name in datasets:
            data = f[ds_name][:].astype(np.float32)

            # Extract wavelength from name (e.g. rhos_561 -> 561 nm)
            wl_str = ds_name.replace(prefix, '')
            try:
                wl = float(wl_str)
            except ValueError:
                wl = 0.0

            # Mask fill values
            data[data == -9999] = np.nan
            data[data == 0] = np.nan  # ACOLITE uses 0 as nodata for reflectance

            arrays.append(data)
            band_names.append(f'{product}_{wl_str} ({wl:.0f} nm)')
            wavelengths.append(wl)

        if not arrays:
            return None

        height, width = arrays[0].shape

        # Get sensor info from attributes
        sensor = attrs.get('sensor', b'unknown')
        if isinstance(sensor, bytes):
            sensor = sensor.decode('utf-8')

        return {
            'arrays': arrays,
            'band_names': band_names,
            'wavelengths': wavelengths,
            'transform': transform,
            'crs': crs,
            'height': height,
            'width': width,
            'sensor': sensor,
            'product': product,
        }


def write_geotiff(data, output_path):
    """Write multiband GeoTIFF."""
    profile = {
        'driver': 'GTiff',
        'dtype': 'float32',
        'width': data['width'],
        'height': data['height'],
        'count': len(data['arrays']),
        'nodata': np.nan,
        'compress': 'deflate',
        'predictor': 2,
        'zlevel': 6,
    }
    if data['crs']:
        profile['crs'] = data['crs']
    if data['transform']:
        profile['transform'] = data['transform']

    with rasterio.open(output_path, 'w', **profile) as dst:
        for i, (arr, name) in enumerate(zip(data['arrays'], data['band_names']), 1):
            dst.write(arr, i)
            dst.set_band_description(i, name)

    print(f"  [GeoTIFF] {os.path.basename(output_path)}  "
          f"({len(data['arrays'])} bands, {data['width']}x{data['height']})")


def write_envi(data, output_base):
    """Write multiband ENVI .dat + .hdr pair."""
    dat_path = output_base + '.dat'
    hdr_path = output_base + '.hdr'

    num_bands = len(data['arrays'])
    height, width = data['height'], data['width']

    # Stack and write binary
    stack = np.zeros((num_bands, height, width), dtype=np.float32)
    for i, arr in enumerate(data['arrays']):
        out = arr.copy()
        out[np.isnan(out)] = -9999.0
        stack[i] = out
    stack.tofile(dat_path)

    # Write header
    with open(hdr_path, 'w') as f:
        f.write('ENVI\n')
        f.write(f'description = {{ACOLITE L2W {data["product"]} - {data["sensor"]}}}\n')
        f.write(f'samples = {width}\n')
        f.write(f'lines = {height}\n')
        f.write(f'bands = {num_bands}\n')
        f.write(f'header offset = 0\n')
        f.write(f'file type = ENVI Standard\n')
        f.write(f'data type = 4\n')  # float32
        f.write(f'interleave = bsq\n')
        f.write(f'byte order = 0\n')
        f.write(f'data ignore value = -9999.0\n')

        if data['transform']:
            t = data['transform']
            f.write(f'map info = {{Arbitrary, 1, 1, {t.c}, {t.f}, '
                    f'{abs(t.a)}, {abs(t.e)}}}\n')

        if data['crs']:
            f.write(f'coordinate system string = {{{data["crs"].to_wkt()}}}\n')

        if data['wavelengths']:
            wl_str = ', '.join(f'{w:.1f}' for w in data['wavelengths'])
            f.write(f'wavelength units = Nanometers\n')
            f.write(f'wavelength = {{\n  {wl_str}}}\n')

        if data['band_names']:
            bn_str = ',\n  '.join(data['band_names'])
            f.write(f'band names = {{\n  {bn_str}}}\n')

    print(f"  [ENVI]    {os.path.basename(dat_path)}  ({num_bands} bands)")


def process_file(filepath, output_dir=None, write_gtiff=True, write_env=True,
                 products=('rhos',)):
    """Process a single ACOLITE L2W file."""
    if output_dir is None:
        output_dir = os.path.join(os.path.dirname(filepath), 'converted')

    basename = os.path.splitext(os.path.basename(filepath))[0]
    print(f"\n{'='*60}")
    print(f"Processing: {os.path.basename(filepath)}")
    print(f"{'='*60}")

    for product in products:
        print(f"\n  Reading {product} bands...")
        data = read_acolite_l2w(filepath, product=product)

        if data is None:
            print(f"  No {product} data found, skipping.")
            continue

        print(f"  Sensor: {data['sensor']}")
        print(f"  Bands: {len(data['arrays'])}")
        print(f"  Wavelengths: {', '.join(f'{w:.0f}' for w in data['wavelengths'])} nm")
        print(f"  Size: {data['width']} x {data['height']}")
        if data['crs']:
            print(f"  CRS: {data['crs']}")

        suffix = f'_{product}' if len(products) > 1 else ''
        out_base = f"{basename}{suffix}"

        if write_gtiff:
            tif_dir = os.path.join(output_dir, 'GeoTIFF')
            os.makedirs(tif_dir, exist_ok=True)
            write_geotiff(data, os.path.join(tif_dir, f"{out_base}.tif"))

        if write_env:
            envi_dir = os.path.join(output_dir, 'ENVI')
            os.makedirs(envi_dir, exist_ok=True)
            write_envi(data, os.path.join(envi_dir, out_base))


def process_directory(input_path, output_dir=None, write_gtiff=True, write_env=True,
                      products=('rhos',)):
    """Process a single file or all .nc files in a directory."""
    if os.path.isfile(input_path):
        process_file(input_path, output_dir, write_gtiff, write_env, products)
    elif os.path.isdir(input_path):
        nc_files = sorted(glob.glob(os.path.join(input_path, '*L2W*.nc')))
        if not nc_files:
            nc_files = sorted(glob.glob(os.path.join(input_path, '*.nc')))
        if not nc_files:
            print("No .nc files found.")
            return
        print(f"Found {len(nc_files)} file(s)")
        for f in nc_files:
            process_file(f, output_dir, write_gtiff, write_env, products)
    else:
        print(f"Error: {input_path} not found")


# ============================================================================
# CLI
# ============================================================================

if __name__ == '__main__':
    if not _HAS_H5PY:
        print("ERROR: h5py is required.  pip install h5py")
        sys.exit(1)
    if not _HAS_RASTERIO:
        print("ERROR: rasterio is required.  conda install -c conda-forge rasterio")
        sys.exit(1)

    if len(sys.argv) < 2:
        print("ACOLITE L2W Converter")
        print("=" * 40)
        print(f"\nUsage: python {os.path.basename(__file__)} <input> [output_dir] [options]")
        print(f"\n  <input> can be a single .nc file or a folder containing .nc files")
        print(f"\nOptions:")
        print(f"  --geotiff    GeoTIFF output only")
        print(f"  --envi       ENVI output only")
        print(f"  --rhot       Extract TOA reflectance (rhot) instead of surface (rhos)")
        print(f"  --both       Extract both rhos and rhot as separate files")
        print(f"\nExamples:")
        print(f"  python {os.path.basename(__file__)} L8_OLI_2020_01_06_L2W.nc")
        print(f"  python {os.path.basename(__file__)} ./acolite_output/ --both")
        sys.exit(0)

    args = sys.argv[1:]
    write_gtiff = True
    write_env = True
    products = ['rhos']
    positional = []

    for a in args:
        al = a.lower()
        if al == '--geotiff':
            write_gtiff = True; write_env = False
        elif al == '--envi':
            write_gtiff = False; write_env = True
        elif al == '--rhot':
            products = ['rhot']
        elif al == '--both':
            products = ['rhos', 'rhot']
        else:
            positional.append(a)

    input_path = positional[0]
    output_dir = positional[1] if len(positional) > 1 else None

    process_directory(input_path, output_dir, write_gtiff, write_env, tuple(products))
    print(f"\nDone.")