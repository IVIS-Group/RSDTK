"""
Landsat Collection 2 Level-1 Converter
========================================
Reads individual band GeoTIFFs + MTL metadata, applies radiometric calibration,
and outputs multiband files. Auto-detects sensor from MTL metadata.

Supported sensors:
  Landsat 8/9 OLI/TIRS (LC08/LC09):
    1. VSWIR (30m): Bands 1-7 + 9 -> TOA Reflectance (8 bands)
    2. TIRS (30m): Bands 10-11 -> At-sensor Radiance (2 bands)
    3. Pan (15m): Band 8 -> TOA Reflectance (1 band)

  Landsat 7 ETM+ (LE07):
    1. VSWIR (30m): Bands 1-5 + 7 -> TOA Reflectance (6 bands)
    2. TIR (60m->30m): Band 6 VCID 1 + 2 -> At-sensor Radiance (2 bands)
    3. Pan (15m): Band 8 -> TOA Reflectance (1 band)

  Landsat 4/5 TM (LT04/LT05):
    1. VSWIR (30m): Bands 1-5 + 7 -> TOA Reflectance (6 bands)
    2. TIR (120m->30m): Band 6 -> At-sensor Radiance (1 band)

  Landsat 1-5 MSS (LM01-LM05):
    1. VSWIR (60m): 4 bands -> TOA Reflectance
       L1-3: Bands 4-7  |  L4-5: Bands 1-4

Dependencies: numpy, rasterio
"""

import os
import re
import sys
import numpy as np

try:
    import rasterio
    from rasterio.crs import CRS
except ImportError:
    print("ERROR: rasterio is required. Install with: pip install rasterio")
    sys.exit(1)


# ============================================================================
# Sensor-specific band configuration
# ============================================================================

SENSOR_CONFIG = {
    'OLI_TIRS': {
        'label': 'Landsat 8/9 OLI/TIRS',
        'vswir_bands': [1, 2, 3, 4, 5, 6, 7, 9],
        'tir_bands': [10, 11],
        'pan_bands': [8],
        'vswir_wavelengths': [0.443, 0.482, 0.561, 0.655, 0.865, 1.609, 2.201, 1.373],
        'tir_wavelengths': [10.895, 12.005],
        'pan_wavelengths': [0.592],
        'vswir_names': [
            'Band 1 - Coastal Aerosol (0.443 um)', 'Band 2 - Blue (0.482 um)',
            'Band 3 - Green (0.561 um)', 'Band 4 - Red (0.655 um)',
            'Band 5 - NIR (0.865 um)', 'Band 6 - SWIR-1 (1.609 um)',
            'Band 7 - SWIR-2 (2.201 um)', 'Band 9 - Cirrus (1.373 um)',
        ],
        'tir_names': ['Band 10 - TIRS-1 (10.895 um)', 'Band 11 - TIRS-2 (12.005 um)'],
        'pan_names': ['Band 8 - Panchromatic (0.592 um)'],
    },
    'ETM': {
        'label': 'Landsat 7 ETM+',
        'vswir_bands': [1, 2, 3, 4, 5, 7],
        'tir_bands': ['6_VCID_1', '6_VCID_2'],
        'pan_bands': [8],
        'vswir_wavelengths': [0.485, 0.560, 0.660, 0.835, 1.650, 2.215],
        'tir_wavelengths': [11.450, 11.450],
        'pan_wavelengths': [0.710],
        'vswir_names': [
            'Band 1 - Blue (0.485 um)', 'Band 2 - Green (0.560 um)',
            'Band 3 - Red (0.660 um)', 'Band 4 - NIR (0.835 um)',
            'Band 5 - SWIR-1 (1.650 um)', 'Band 7 - SWIR-2 (2.215 um)',
        ],
        'tir_names': [
            'Band 6 VCID 1 - TIR Low Gain (11.450 um)',
            'Band 6 VCID 2 - TIR High Gain (11.450 um)',
        ],
        'pan_names': ['Band 8 - Panchromatic (0.710 um)'],
    },
    'TM': {
        'label': 'Landsat 4/5 TM',
        'vswir_bands': [1, 2, 3, 4, 5, 7],
        'tir_bands': [6],
        'pan_bands': [],
        'vswir_wavelengths': [0.485, 0.560, 0.660, 0.830, 1.650, 2.215],
        'tir_wavelengths': [11.450],
        'pan_wavelengths': [],
        'vswir_names': [
            'Band 1 - Blue (0.485 um)', 'Band 2 - Green (0.560 um)',
            'Band 3 - Red (0.660 um)', 'Band 4 - NIR (0.830 um)',
            'Band 5 - SWIR-1 (1.650 um)', 'Band 7 - SWIR-2 (2.215 um)',
        ],
        'tir_names': ['Band 6 - TIR (11.450 um)'],
        'pan_names': [],
    },
    'MSS_EARLY': {
        'label': 'Landsat 1-3 MSS',
        'vswir_bands': [4, 5, 6, 7],
        'tir_bands': [],
        'pan_bands': [],
        'vswir_wavelengths': [0.550, 0.650, 0.750, 0.950],
        'tir_wavelengths': [],
        'pan_wavelengths': [],
        'vswir_names': [
            'Band 4 - Green (0.550 um)', 'Band 5 - Red (0.650 um)',
            'Band 6 - NIR-1 (0.750 um)', 'Band 7 - NIR-2 (0.950 um)',
        ],
        'tir_names': [],
        'pan_names': [],
    },
    'MSS_LATE': {
        'label': 'Landsat 4/5 MSS',
        'vswir_bands': [1, 2, 3, 4],
        'tir_bands': [],
        'pan_bands': [],
        'vswir_wavelengths': [0.550, 0.650, 0.750, 0.950],
        'tir_wavelengths': [],
        'pan_wavelengths': [],
        'vswir_names': [
            'Band 1 - Green (0.550 um)', 'Band 2 - Red (0.650 um)',
            'Band 3 - NIR-1 (0.750 um)', 'Band 4 - NIR-2 (0.950 um)',
        ],
        'tir_names': [],
        'pan_names': [],
    },
}


def detect_sensor(metadata):
    """Detect the Landsat sensor type from parsed MTL metadata."""
    root = metadata.get('LANDSAT_METADATA_FILE', metadata)
    img_attr = root.get('IMAGE_ATTRIBUTES', {})

    spacecraft = img_attr.get('SPACECRAFT_ID', '')
    sensor_id = img_attr.get('SENSOR_ID', '')

    if sensor_id == 'MSS':
        sat_num = 0
        match = re.search(r'LANDSAT_(\d+)', spacecraft)
        if match:
            sat_num = int(match.group(1))
        return 'MSS_EARLY' if sat_num <= 3 else 'MSS_LATE'
    elif sensor_id == 'TM':
        return 'TM'
    elif sensor_id in ('ETM', 'ETM+'):
        return 'ETM'
    elif sensor_id in ('OLI_TIRS', 'OLI', 'TIRS'):
        return 'OLI_TIRS'

    # Fallback: detect from product ID
    prod_contents = root.get('PRODUCT_CONTENTS', {})
    product_id = prod_contents.get('LANDSAT_PRODUCT_ID', '')
    if product_id.startswith('LC0'):
        return 'OLI_TIRS'
    elif product_id.startswith('LE07'):
        return 'ETM'
    elif product_id.startswith('LT0'):
        return 'TM'
    elif product_id.startswith('LM0'):
        match = re.match(r'LM0(\d)', product_id)
        if match and int(match.group(1)) <= 3:
            return 'MSS_EARLY'
        return 'MSS_LATE'

    return 'OLI_TIRS'


# ============================================================================
# MTL metadata parser
# ============================================================================

def parse_mtl(mtl_path):
    """Parse a Landsat MTL.txt file into a nested dictionary."""
    metadata = {}
    group_stack = [metadata]

    with open(mtl_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line == 'END':
                continue
            if line.startswith('GROUP'):
                _, group_name = line.split('=')
                group_name = group_name.strip()
                new_group = {}
                group_stack[-1][group_name] = new_group
                group_stack.append(new_group)
            elif line.startswith('END_GROUP'):
                group_stack.pop()
            elif '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip().strip('"')
                try:
                    if '.' in value or 'E' in value.upper():
                        value = float(value)
                    else:
                        value = int(value)
                except (ValueError, TypeError):
                    pass
                group_stack[-1][key] = value

    return metadata


def extract_calibration(metadata):
    """
    Extract all calibration coefficients from parsed MTL metadata.
    Handles standard band numbering and ETM+ VCID naming.
    """
    cal = {
        'sun_elevation': None,
        'radiance_mult': {}, 'radiance_add': {},
        'reflectance_mult': {}, 'reflectance_add': {},
        'k1': {}, 'k2': {},
    }

    root = metadata.get('LANDSAT_METADATA_FILE', metadata)
    img_attr = root.get('IMAGE_ATTRIBUTES', {})
    cal['sun_elevation'] = img_attr.get('SUN_ELEVATION')

    rescaling = root.get('LEVEL1_RADIOMETRIC_RESCALING', {})
    for key, value in rescaling.items():
        if not isinstance(value, (int, float)):
            continue

        # Standard bands: RADIANCE_MULT_BAND_1
        m = re.match(r'RADIANCE_MULT_BAND_(\d+)$', key)
        if m:
            cal['radiance_mult'][int(m.group(1))] = value
            continue
        m = re.match(r'RADIANCE_ADD_BAND_(\d+)$', key)
        if m:
            cal['radiance_add'][int(m.group(1))] = value
            continue
        m = re.match(r'REFLECTANCE_MULT_BAND_(\d+)$', key)
        if m:
            cal['reflectance_mult'][int(m.group(1))] = value
            continue
        m = re.match(r'REFLECTANCE_ADD_BAND_(\d+)$', key)
        if m:
            cal['reflectance_add'][int(m.group(1))] = value
            continue

        # ETM+ VCID: RADIANCE_MULT_BAND_6_VCID_1
        m = re.match(r'RADIANCE_MULT_BAND_(\d+_VCID_\d+)', key)
        if m:
            cal['radiance_mult'][m.group(1)] = value
            continue
        m = re.match(r'RADIANCE_ADD_BAND_(\d+_VCID_\d+)', key)
        if m:
            cal['radiance_add'][m.group(1)] = value
            continue

    thermal = root.get('LEVEL1_THERMAL_CONSTANTS', {})
    for key, value in thermal.items():
        if not isinstance(value, (int, float)):
            continue
        m = re.match(r'K1_CONSTANT_BAND_(\S+)', key)
        if m:
            bk = m.group(1)
            try:
                cal['k1'][int(bk)] = value
            except ValueError:
                cal['k1'][bk] = value
            continue
        m = re.match(r'K2_CONSTANT_BAND_(\S+)', key)
        if m:
            bk = m.group(1)
            try:
                cal['k2'][int(bk)] = value
            except ValueError:
                cal['k2'][bk] = value

    return cal


# ============================================================================
# Scene discovery
# ============================================================================

def discover_scenes(input_dir):
    """Scan for Landsat L1 scenes (identified by _MTL.txt).
    When multiple MTL files exist in one directory, prefer L1TP/L1GT."""
    scenes = []

    def pick_l1_mtl(mtl_list):
        """From a list of MTL filenames, prefer L1TP/L1GT over L2."""
        l1 = [f for f in mtl_list if '_L1TP_' in f or '_L1GT_' in f]
        if l1:
            return l1[0]
        return mtl_list[0]

    mtl_files = [f for f in os.listdir(input_dir) if f.endswith('_MTL.txt')]
    if mtl_files:
        chosen = pick_l1_mtl(mtl_files)
        mtl_path = os.path.join(input_dir, chosen)
        scene_id = chosen.replace('_MTL.txt', '')
        scenes.append((input_dir, scene_id, mtl_path))
        return scenes

    for entry in sorted(os.listdir(input_dir)):
        subdir = os.path.join(input_dir, entry)
        if os.path.isdir(subdir):
            mtl_files = [f for f in os.listdir(subdir)
                         if f.endswith('_MTL.txt')]
            if mtl_files:
                chosen = pick_l1_mtl(mtl_files)
                mtl_path = os.path.join(subdir, chosen)
                scene_id = chosen.replace('_MTL.txt', '')
                scenes.append((subdir, scene_id, mtl_path))

    return scenes


def find_band_file(scene_dir, scene_id, band_id):
    """
    Find the band file for a given band identifier.
    band_id can be int (1) or string ('6_VCID_1').
    """
    suffix = f'_B{band_id}.TIF'
    for f in os.listdir(scene_dir):
        if f.upper().endswith(suffix.upper()):
            base = f.upper().replace('.TIF', '')
            if base.endswith(f'_B{str(band_id).upper()}'):
                return os.path.join(scene_dir, f)
    return None


# ============================================================================
# Band reading and calibration
# ============================================================================

def read_and_calibrate_reflectance(filepath, mult, add, sun_elevation,
                                   fill_value=0):
    """TOA Reflectance = (MULT * DN + ADD) / sin(SUN_ELEVATION)"""
    with rasterio.open(filepath) as src:
        dn = src.read(1)
        profile = src.profile.copy()
    data = np.float32(dn)
    fill_mask = (dn == fill_value)
    sun_rad = np.sin(np.radians(sun_elevation))
    data = (mult * data + add) / sun_rad
    data[fill_mask] = np.nan
    return data, profile


def read_and_calibrate_radiance(filepath, mult, add, fill_value=0):
    """Radiance = MULT * DN + ADD"""
    with rasterio.open(filepath) as src:
        dn = src.read(1)
        profile = src.profile.copy()
    data = np.float32(dn)
    fill_mask = (dn == fill_value)
    data = mult * data + add
    data[fill_mask] = np.nan
    return data, profile


# ============================================================================
# ENVI header writing
# ============================================================================

def build_envi_header(envi_path, height, width, num_bands, dtype, crs,
                      transform, wavelengths=None, band_names=None,
                      description='', interleave='bsq'):
    """Write a custom ENVI header file with wavelength metadata."""
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
                f.write(f'map info = {{Geographic Lat/Lon, 1, 1, {ul_x}, '
                        f'{ul_y}, {x_size}, {y_size}, WGS-84}}\n')
            elif crs and crs.is_projected:
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
            f.write('wavelength units = Micrometers\n')
            wl_str = ', '.join(f'{w:.3f}' for w in wavelengths)
            f.write(f'wavelength = {{\n  {wl_str}}}\n')

        if band_names is not None:
            bn_str = ',\n  '.join(band_names)
            f.write(f'band names = {{\n  {bn_str}}}\n')

        f.write('data ignore value = nan\n')


# ============================================================================
# Output writing
# ============================================================================

def write_multiband(band_data_list, output_path, ref_profile, band_names,
                    wavelengths=None, description='',
                    write_geotiff=True, write_envi=True):
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
        print(f"    Writing GeoTIFF: {os.path.basename(gtiff_path)}")
        with rasterio.open(gtiff_path, 'w', **profile) as dst:
            for i in range(num_bands):
                dst.write(stack[i], i + 1)
                dst.set_band_description(i + 1, band_names[i])
        output_files.append(gtiff_path)

    if write_envi:
        envi_dat = output_path + '.dat'
        envi_hdr = output_path + '.hdr'
        print(f"    Writing ENVI: {os.path.basename(envi_dat)}")
        stack.tofile(envi_dat)
        build_envi_header(
            envi_hdr, height=height, width=width, num_bands=num_bands,
            dtype=np.float32, crs=crs, transform=transform,
            wavelengths=wavelengths, band_names=band_names,
            description=description, interleave='bsq'
        )
        output_files.append(envi_dat)
        output_files.append(envi_hdr)

    return output_files


# ============================================================================
# Main processing
# ============================================================================

def process_scene(scene_dir, scene_id, mtl_path, output_dir,
                  write_geotiff=True, write_envi=True, overwrite=True):
    """Process a single Landsat L1 scene. Auto-detects sensor."""
    all_outputs = []

    print(f"  Parsing MTL: {os.path.basename(mtl_path)}")
    metadata = parse_mtl(mtl_path)
    cal = extract_calibration(metadata)

    sensor_key = detect_sensor(metadata)
    config = SENSOR_CONFIG[sensor_key]
    print(f"  Sensor: {config['label']}")

    sun_elev = cal['sun_elevation']
    if sun_elev is None:
        print("  ERROR: SUN_ELEVATION not found in MTL. Skipping scene.")
        return all_outputs
    print(f"  Sun elevation: {sun_elev:.4f} degrees")

    os.makedirs(output_dir, exist_ok=True)

    # --- VSWIR: TOA Reflectance ---
    if config['vswir_bands']:
        print(f"\n  === Processing VSWIR (TOA Reflectance) ===")
        vswir_data = []
        ref_profile = None
        actual_names = []
        actual_wl = []

        for i, band_id in enumerate(config['vswir_bands']):
            filepath = find_band_file(scene_dir, scene_id, band_id)
            if filepath is None:
                print(f"    WARNING: Band {band_id} not found, skipping.")
                continue

            mult = cal['reflectance_mult'].get(band_id)
            add = cal['reflectance_add'].get(band_id)

            if mult is None or add is None:
                # Fallback to radiance for bands without reflectance coeffs
                rmult = cal['radiance_mult'].get(band_id)
                radd = cal['radiance_add'].get(band_id)
                if rmult is not None and radd is not None:
                    print(f"    Reading Band {band_id} as Radiance "
                          f"(MULT={rmult}, ADD={radd})...")
                    data, profile = read_and_calibrate_radiance(
                        filepath, rmult, radd)
                else:
                    print(f"    WARNING: No calibration for Band {band_id}.")
                    continue
            else:
                print(f"    Reading Band {band_id} (MULT={mult}, ADD={add})...")
                data, profile = read_and_calibrate_reflectance(
                    filepath, mult, add, sun_elev)

            vswir_data.append(data)
            actual_names.append(config['vswir_names'][i]
                                if i < len(config['vswir_names'])
                                else f'Band {band_id}')
            actual_wl.append(config['vswir_wavelengths'][i]
                             if i < len(config['vswir_wavelengths'])
                             else 0.0)
            if ref_profile is None:
                ref_profile = profile
                print(f"    Image: {profile['width']} x {profile['height']}")

        if vswir_data and ref_profile:
            vswir_base = os.path.join(output_dir, f"{scene_id}_VSWIR_TOARef")
            if not overwrite and os.path.isfile(vswir_base + '.tif'):
                print("    Output exists, skipping (overwrite=False)")
            else:
                outputs = write_multiband(
                    vswir_data, vswir_base, ref_profile,
                    band_names=actual_names, wavelengths=actual_wl,
                    description=f'{config["label"]} L1 TOA Reflectance - '
                                f'VSWIR ({scene_id})',
                    write_geotiff=write_geotiff, write_envi=write_envi)
                all_outputs.extend(outputs)

    # --- TIR: At-sensor Radiance ---
    if config['tir_bands']:
        print(f"\n  === Processing TIR (At-sensor Radiance) ===")
        tir_data = []
        tir_profile = None
        actual_tir_names = []
        actual_tir_wl = []

        for i, band_id in enumerate(config['tir_bands']):
            filepath = find_band_file(scene_dir, scene_id, band_id)
            if filepath is None:
                print(f"    WARNING: Band {band_id} not found, skipping.")
                continue

            mult = cal['radiance_mult'].get(band_id)
            add = cal['radiance_add'].get(band_id)
            if mult is None or add is None:
                print(f"    WARNING: No radiance coefficients for Band "
                      f"{band_id}, skipping.")
                continue

            print(f"    Reading Band {band_id} (MULT={mult}, ADD={add})...")
            data, profile = read_and_calibrate_radiance(filepath, mult, add)
            tir_data.append(data)
            actual_tir_names.append(config['tir_names'][i]
                                    if i < len(config['tir_names'])
                                    else f'Band {band_id}')
            actual_tir_wl.append(config['tir_wavelengths'][i]
                                 if i < len(config['tir_wavelengths'])
                                 else 0.0)
            if tir_profile is None:
                tir_profile = profile

        if tir_data and tir_profile:
            tir_base = os.path.join(output_dir, f"{scene_id}_TIR_Radiance")
            if not overwrite and os.path.isfile(tir_base + '.tif'):
                print("    Output exists, skipping (overwrite=False)")
            else:
                outputs = write_multiband(
                    tir_data, tir_base, tir_profile,
                    band_names=actual_tir_names, wavelengths=actual_tir_wl,
                    description=f'{config["label"]} L1 At-sensor Radiance - '
                                f'TIR ({scene_id})',
                    write_geotiff=write_geotiff, write_envi=write_envi)
                all_outputs.extend(outputs)

    # --- Pan: TOA Reflectance ---
    if config['pan_bands']:
        print(f"\n  === Processing Panchromatic (TOA Reflectance) ===")
        pan_data = []
        pan_profile = None

        for band_id in config['pan_bands']:
            filepath = find_band_file(scene_dir, scene_id, band_id)
            if filepath is None:
                print(f"    WARNING: Band {band_id} not found, skipping.")
                continue
            mult = cal['reflectance_mult'].get(band_id)
            add = cal['reflectance_add'].get(band_id)
            if mult is None or add is None:
                print(f"    WARNING: No reflectance coefficients for Band "
                      f"{band_id}, skipping.")
                continue
            print(f"    Reading Band {band_id} (MULT={mult}, ADD={add})...")
            data, profile = read_and_calibrate_reflectance(
                filepath, mult, add, sun_elev)
            pan_data.append(data)
            if pan_profile is None:
                pan_profile = profile
                print(f"    Pan: {profile['width']} x {profile['height']}")

        if pan_data and pan_profile:
            pan_base = os.path.join(output_dir, f"{scene_id}_Pan_TOARef")
            if not overwrite and os.path.isfile(pan_base + '.tif'):
                print("    Output exists, skipping (overwrite=False)")
            else:
                outputs = write_multiband(
                    pan_data, pan_base, pan_profile,
                    band_names=config['pan_names'][:len(pan_data)],
                    wavelengths=config['pan_wavelengths'][:len(pan_data)],
                    description=f'{config["label"]} L1 TOA Reflectance - '
                                f'Panchromatic ({scene_id})',
                    write_geotiff=write_geotiff, write_envi=write_envi)
                all_outputs.extend(outputs)

    return all_outputs


def process_directory(input_dir, output_dir=None, write_geotiff=True,
                      write_envi=True, overwrite=True):
    """Process all Landsat L1 scenes found in the input directory."""
    if output_dir is None:
        output_dir = os.path.join(input_dir, 'converted')

    print(f"Scanning for Landsat L1 scenes in: {input_dir}")
    scenes = discover_scenes(input_dir)

    if not scenes:
        print("No Landsat L1 scenes found.")
        return

    print(f"Found {len(scenes)} scene(s):\n")

    all_outputs = []
    for i, (scene_dir, scene_id, mtl_path) in enumerate(scenes, 1):
        print(f"{'='*70}")
        print(f"Scene {i}/{len(scenes)}: {scene_id}")
        print(f"{'='*70}")
        try:
            outputs = process_scene(
                scene_dir, scene_id, mtl_path, output_dir,
                write_geotiff=write_geotiff, write_envi=write_envi,
                overwrite=overwrite)
            all_outputs.extend(outputs)
        except Exception as e:
            print(f"  ERROR processing scene {scene_id}: {e}")
            import traceback
            traceback.print_exc()
        print()

    print(f"Processing complete. {len(all_outputs)} output file(s) in: "
          f"{output_dir}")


# ============================================================================
# CLI
# ============================================================================

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Landsat Collection 2 Level-1 Converter")
        print("=" * 50)
        print(f"\nUsage: python {os.path.basename(__file__)} <input_dir> "
              f"[output_dir] [options]")
        print(f"\nSupported sensors (auto-detected from MTL):")
        print(f"  Landsat 8/9 OLI/TIRS  |  Landsat 7 ETM+")
        print(f"  Landsat 4/5 TM        |  Landsat 1-5 MSS")
        print(f"\nOptions: --geotiff, --envi (default: both)")
        sys.exit(0)

    args = sys.argv[1:]
    write_geotiff = True
    write_envi = True
    positional = []

    for arg in args:
        if arg.lower() == '--geotiff':
            write_geotiff, write_envi = True, False
        elif arg.lower() == '--envi':
            write_geotiff, write_envi = False, True
        else:
            positional.append(arg)

    input_dir = positional[0]
    output_dir = positional[1] if len(positional) > 1 else None

    if not os.path.isdir(input_dir):
        print(f"Error: {input_dir} does not exist")
        sys.exit(1)

    process_directory(input_dir, output_dir,
                      write_geotiff=write_geotiff, write_envi=write_envi)
