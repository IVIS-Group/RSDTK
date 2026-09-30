"""
Sentinel-2 L2A Surface Reflectance Converter
Reads Sentinel-2 L2A .SAFE folders, applies BOA reflectance scale factors,
and outputs multiband GeoTIFF and/or ENVI files:
  1. R10m (4 bands): B02, B03, B04, B08 at 10m
  2. R20m (10 bands): B01-B07, B8A, B11, B12 at 20m
  3. SCL (1 band): Scene Classification Layer at 20m

Supports Sentinel-2A and Sentinel-2B.
Scale: Reflectance = (DN + BOA_ADD_OFFSET) / BOA_QUANTIFICATION_VALUE

Band Wavelengths (nm):
  B01: 443 (Coastal)    B02: 490 (Blue)      B03: 560 (Green)
  B04: 665 (Red)        B05: 705 (Red Edge 1) B06: 740 (Red Edge 2)
  B07: 783 (Red Edge 3) B08: 842 (NIR)       B8A: 865 (NIR narrow)
  B09: 945 (Water Vap)  B11: 1610 (SWIR-1)   B12: 2190 (SWIR-2)

Dependencies: rasterio, numpy
"""

import os
import re
import sys
import xml.etree.ElementTree as ET
import numpy as np

try:
    import rasterio
    from rasterio.crs import CRS
except ImportError:
    print("ERROR: rasterio is required. Install with: pip install rasterio")
    sys.exit(1)


# ============================================================================
# Constants
# ============================================================================

# Band definitions: (band_id, center_wavelength_nm, name)
R10M_BANDS = [
    ('B02', 490, 'B02 - Blue (490 nm)'),
    ('B03', 560, 'B03 - Green (560 nm)'),
    ('B04', 665, 'B04 - Red (665 nm)'),
    ('B08', 842, 'B08 - NIR (842 nm)'),
]

R20M_BANDS = [
    ('B01', 443, 'B01 - Coastal Aerosol (443 nm)'),
    ('B02', 490, 'B02 - Blue (490 nm)'),
    ('B03', 560, 'B03 - Green (560 nm)'),
    ('B04', 665, 'B04 - Red (665 nm)'),
    ('B05', 705, 'B05 - Red Edge 1 (705 nm)'),
    ('B06', 740, 'B06 - Red Edge 2 (740 nm)'),
    ('B07', 783, 'B07 - Red Edge 3 (783 nm)'),
    ('B8A', 865, 'B8A - NIR narrow (865 nm)'),
    ('B11', 1610, 'B11 - SWIR-1 (1610 nm)'),
    ('B12', 2190, 'B12 - SWIR-2 (2190 nm)'),
]

# Default scale factors (can be overridden from metadata)
DEFAULT_BOA_QUANTIFICATION = 10000
DEFAULT_BOA_OFFSET = -1000


# ============================================================================
# Scene discovery
# ============================================================================

def discover_scenes(input_dir):
    """
    Find Sentinel-2 L2A .SAFE folders.
    Checks input_dir itself and immediate subdirectories.
    Returns list of (safe_path, scene_id) tuples.
    """
    scenes = []

    def check_safe(path):
        if path.endswith('.SAFE') and os.path.isdir(path):
            mtd = os.path.join(path, 'MTD_MSIL2A.xml')
            if os.path.exists(mtd):
                scene_id = os.path.basename(path).replace('.SAFE', '')
                scenes.append((path, scene_id))

    # Check input_dir itself
    check_safe(input_dir)

    # Check subdirectories
    for entry in os.listdir(input_dir):
        check_safe(os.path.join(input_dir, entry))

    return scenes


# ============================================================================
# Metadata parsing
# ============================================================================

def parse_metadata(safe_path):
    """
    Parse MTD_MSIL2A.xml to extract scale factors.
    Returns dict with quantification value and offset.
    """
    mtd_path = os.path.join(safe_path, 'MTD_MSIL2A.xml')

    cal = {
        'boa_quantification': DEFAULT_BOA_QUANTIFICATION,
        'boa_offset': DEFAULT_BOA_OFFSET,
    }

    try:
        tree = ET.parse(mtd_path)
        root = tree.getroot()

        # Handle namespace
        ns = ''
        if root.tag.startswith('{'):
            ns = root.tag.split('}')[0] + '}'

        # Find quantification value
        for elem in root.iter():
            tag = elem.tag.replace(ns, '')
            if tag == 'BOA_QUANTIFICATION_VALUE' and elem.text:
                cal['boa_quantification'] = float(elem.text)
            elif tag == 'BOA_ADD_OFFSET' and elem.text:
                cal['boa_offset'] = int(elem.text)
                break  # All offsets are the same

    except Exception as e:
        print(f"  WARNING: Could not parse metadata: {e}")
        print(f"  Using defaults: quantification={cal['boa_quantification']}, "
              f"offset={cal['boa_offset']}")

    print(f"  BOA Quantification: {cal['boa_quantification']}")
    print(f"  BOA Offset: {cal['boa_offset']}")

    return cal


def find_granule_path(safe_path):
    """Find the granule IMG_DATA path within the .SAFE folder."""
    granule_dir = os.path.join(safe_path, 'GRANULE')
    if not os.path.isdir(granule_dir):
        raise FileNotFoundError(f"GRANULE directory not found in {safe_path}")

    for entry in os.listdir(granule_dir):
        img_path = os.path.join(granule_dir, entry, 'IMG_DATA')
        if os.path.isdir(img_path):
            return img_path

    raise FileNotFoundError(f"No IMG_DATA directory found in {granule_dir}")


def find_band_file(img_data_path, resolution_dir, band_id):
    """
    Find a specific band JP2 file in the resolution folder.
    Handles both _B02_10m.jp2 and _B8A_20m.jp2 naming.
    """
    res_path = os.path.join(img_data_path, resolution_dir)
    if not os.path.isdir(res_path):
        return None

    # Match pattern: *_B02_10m.jp2 or *_B8A_20m.jp2
    for f in os.listdir(res_path):
        if f.lower().endswith('.jp2'):
            # Extract band ID from filename
            parts = f.replace('.jp2', '').split('_')
            for part in parts:
                if part == band_id:
                    return os.path.join(res_path, f)

    return None


def find_scl_file(img_data_path):
    """Find the Scene Classification Layer file."""
    for res in ['R20m', 'R60m']:
        res_path = os.path.join(img_data_path, res)
        if os.path.isdir(res_path):
            for f in os.listdir(res_path):
                if 'SCL' in f and f.lower().endswith('.jp2'):
                    return os.path.join(res_path, f)
    return None


# ============================================================================
# Band reading and calibration
# ============================================================================

def read_and_calibrate(filepath, quantification, offset, fill_value=0):
    """
    Read a JP2 band and convert to surface reflectance.
    Reflectance = (DN + offset) / quantification
    """
    with rasterio.open(filepath) as src:
        dn = src.read(1)
        profile = src.profile.copy()

    data = dn.astype(np.float32)
    fill_mask = (dn == fill_value)

    data = (data + offset) / quantification

    data[fill_mask] = np.nan

    return data, profile


def read_scl(filepath):
    """Read Scene Classification Layer (no scaling needed)."""
    with rasterio.open(filepath) as src:
        data = src.read(1).astype(np.float32)
        profile = src.profile.copy()

    # SCL values: 0=nodata, 1=saturated, 2=dark, 3=shadow,
    # 4=vegetation, 5=bare soil, 6=water, 7=cloud low prob,
    # 8=cloud med prob, 9=cloud high prob, 10=thin cirrus,
    # 11=snow/ice
    return data, profile


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

            if crs and crs.is_projected:
                epsg = crs.to_epsg()
                if epsg and 32600 < epsg < 32661:
                    zone = epsg - 32600
                    f.write(f'map info = {{UTM, 1, 1, {ul_x}, {ul_y}, '
                            f'{x_size}, {y_size}, {zone}, North, WGS-84}}\n')
                elif epsg and 32700 < epsg < 32761:
                    zone = epsg - 32700
                    f.write(f'map info = {{UTM, 1, 1, {ul_x}, {ul_y}, '
                            f'{x_size}, {y_size}, {zone}, South, WGS-84}}\n')

        if crs is not None:
            f.write(f'coordinate system string = {{{crs.to_wkt()}}}\n')

        if wavelengths is not None:
            f.write('wavelength units = Nanometers\n')
            wl_str = ', '.join(f'{w:.1f}' for w in wavelengths)
            f.write(f'wavelength = {{\n  {wl_str}}}\n')

        if band_names is not None:
            bn_str = ',\n  '.join(band_names)
            f.write(f'band names = {{\n  {bn_str}}}\n')

    return hdr_path


# ============================================================================
# Output writing
# ============================================================================

def write_multiband(band_data_list, output_path, ref_profile, band_names,
                    wavelengths=None, description='',
                    write_geotiff=True, write_envi=True, overwrite=True):
    """Write list of 2D arrays as multiband GeoTIFF and/or ENVI."""
    output_files = []
    num_bands = len(band_data_list)
    height, width = band_data_list[0].shape

    stack = np.zeros((num_bands, height, width), dtype=np.float32)
    for i, data in enumerate(band_data_list):
        stack[i] = data

    crs = ref_profile.get('crs')
    transform = ref_profile.get('transform')

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

def process_scene(safe_path, scene_id, output_dir,
                  write_geotiff=True, write_envi=True, overwrite=True):
    """Process a single Sentinel-2 L2A scene."""
    all_outputs = []

    # Parse metadata
    print(f"  Parsing metadata...")
    cal = parse_metadata(safe_path)

    # Find granule path
    img_data = find_granule_path(safe_path)
    print(f"  IMG_DATA: {img_data}")

    os.makedirs(output_dir, exist_ok=True)

    quant = cal['boa_quantification']
    offset = cal['boa_offset']

    # --- R10m bands ---
    print(f"\n  === Processing R10m Bands ===")
    r10_data = []
    r10_profile = None
    r10_band_names = []
    r10_wavelengths = []

    for band_id, wl, name in R10M_BANDS:
        filepath = find_band_file(img_data, 'R10m', band_id)
        if filepath is None:
            print(f"    WARNING: {band_id} not found, skipping")
            continue

        print(f"    Reading {band_id} ({wl} nm)...")
        data, profile = read_and_calibrate(filepath, quant, offset)
        r10_data.append(data)
        r10_band_names.append(name)
        r10_wavelengths.append(wl)

        if r10_profile is None:
            r10_profile = profile
            print(f"    Dimensions: {profile['width']} x {profile['height']} @ 10m")

    if r10_data and r10_profile:
        r10_base = os.path.join(output_dir, f"{scene_id}_R10m_SR")
        outputs = write_multiband(
            r10_data, r10_base, r10_profile,
            band_names=r10_band_names,
            wavelengths=r10_wavelengths,
            description=f'Sentinel-2 L2A Surface Reflectance R10m ({scene_id})',
            write_geotiff=write_geotiff, write_envi=write_envi,
            overwrite=overwrite
        )
        all_outputs.extend(outputs)

    # --- R20m bands ---
    print(f"\n  === Processing R20m Bands ===")
    r20_data = []
    r20_profile = None
    r20_band_names = []
    r20_wavelengths = []

    for band_id, wl, name in R20M_BANDS:
        filepath = find_band_file(img_data, 'R20m', band_id)
        if filepath is None:
            print(f"    WARNING: {band_id} not found at R20m, skipping")
            continue

        print(f"    Reading {band_id} ({wl} nm)...")
        data, profile = read_and_calibrate(filepath, quant, offset)
        r20_data.append(data)
        r20_band_names.append(name)
        r20_wavelengths.append(wl)

        if r20_profile is None:
            r20_profile = profile
            print(f"    Dimensions: {profile['width']} x {profile['height']} @ 20m")

    if r20_data and r20_profile:
        r20_base = os.path.join(output_dir, f"{scene_id}_R20m_SR")
        outputs = write_multiband(
            r20_data, r20_base, r20_profile,
            band_names=r20_band_names,
            wavelengths=r20_wavelengths,
            description=f'Sentinel-2 L2A Surface Reflectance R20m ({scene_id})',
            write_geotiff=write_geotiff, write_envi=write_envi,
            overwrite=overwrite
        )
        all_outputs.extend(outputs)

    # --- SCL ---
    print(f"\n  === Processing Scene Classification ===")
    scl_path = find_scl_file(img_data)

    if scl_path:
        print(f"    Reading SCL...")
        scl_data, scl_profile = read_scl(scl_path)
        print(f"    Dimensions: {scl_profile['width']} x {scl_profile['height']}")

        scl_base = os.path.join(output_dir, f"{scene_id}_SCL")
        outputs = write_multiband(
            [scl_data], scl_base, scl_profile,
            band_names=['Scene Classification Layer'],
            description=f'Sentinel-2 L2A Scene Classification ({scene_id})',
            write_geotiff=write_geotiff, write_envi=write_envi,
            overwrite=overwrite
        )
        all_outputs.extend(outputs)
    else:
        print(f"    WARNING: SCL file not found")

    return all_outputs


def process_directory(input_dir, output_dir=None,
                      write_geotiff=True, write_envi=True,
                      overwrite=True):
    """Process all Sentinel-2 L2A scenes in a directory."""
    if output_dir is None:
        output_dir = os.path.join(input_dir, 'converted')

    print(f"Scanning for Sentinel-2 L2A scenes in: {input_dir}")
    scenes = discover_scenes(input_dir)

    if not scenes:
        print("No Sentinel-2 L2A .SAFE folders found.")
        return

    print(f"Found {len(scenes)} scene(s):\n")

    all_outputs = []
    for i, (safe_path, scene_id) in enumerate(scenes, 1):
        print(f"{'='*70}")
        print(f"Scene {i}/{len(scenes)}: {scene_id}")
        print(f"{'='*70}")

        try:
            outputs = process_scene(
                safe_path, scene_id, output_dir,
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
        print("Sentinel-2 L2A Surface Reflectance Converter")
        print("=" * 50)
        print(f"\nUsage: python {os.path.basename(__file__)} <input_dir> [output_dir] [options]")
        print(f"\nOutput format:")
        print(f"  --geotiff    GeoTIFF only")
        print(f"  --envi       ENVI only")
        print(f"  (default: both)")
        print(f"\nThe input can be a .SAFE folder directly or a directory")
        print(f"containing one or more .SAFE folders.")
        print(f"\nOutputs per scene:")
        print(f"  1. R10m Surface Reflectance (4 bands: B02, B03, B04, B08)")
        print(f"  2. R20m Surface Reflectance (10 bands: B01-B07, B8A, B11, B12)")
        print(f"  3. Scene Classification Layer (1 band)")
        print(f"\nSupports Sentinel-2A and Sentinel-2B (S2A/S2B).")
        sys.exit(0)

    args = sys.argv[1:]
    write_geotiff = True
    write_envi = True
    positional = []

    for arg in args:
        a = arg.lower()
        if a == '--geotiff':
            write_geotiff = True; write_envi = False
        elif a == '--envi':
            write_geotiff = False; write_envi = True
        else:
            positional.append(arg)

    input_dir = positional[0]
    output_dir = positional[1] if len(positional) > 1 else None

    if not os.path.isdir(input_dir):
        print(f"Error: Input directory does not exist: {input_dir}")
        sys.exit(1)

    print(f"Output formats: {'GeoTIFF' if write_geotiff else ''}"
          f"{' + ' if write_geotiff and write_envi else ''}"
          f"{'ENVI' if write_envi else ''}")
    print()

    process_directory(input_dir, output_dir,
                      write_geotiff=write_geotiff, write_envi=write_envi)