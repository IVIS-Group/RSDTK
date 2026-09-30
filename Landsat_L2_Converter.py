"""
Landsat Collection 2 Level-2 (L2SP) Converter
================================================
Reads individual band GeoTIFFs + MTL metadata, applies scale factors,
and outputs multiband files. Auto-detects sensor from MTL metadata.

Supported sensors:
  Landsat 8/9 OLI/TIRS (LC08/LC09):
    1. VSWIR (30m): SR_B1-7 -> Surface Reflectance (7 bands)
    2. TIRS (30m): ST_B10 -> Surface Temperature in Kelvin (1 band)
    3. ST Ancillary (30m): TRAD, URAD, DRAD, ATRAN, EMIS, EMSD, CDIST

  Landsat 7 ETM+ (LE07):
    1. VSWIR (30m): SR_B1-5 + SR_B7 -> Surface Reflectance (6 bands)
    2. TIR (30m): ST_B6 -> Surface Temperature in Kelvin (1 band)
    3. ST Ancillary (30m): TRAD, URAD, DRAD, ATRAN, EMIS, EMSD, CDIST

  Landsat 4/5 TM (LT04/LT05):
    1. VSWIR (30m): SR_B1-5 + SR_B7 -> Surface Reflectance (6 bands)
    2. TIR (30m): ST_B6 -> Surface Temperature in Kelvin (1 band)
    3. ST Ancillary (30m): TRAD, URAD, DRAD, ATRAN, EMIS, EMSD, CDIST

Note: MSS (LM01-LM05) has no Level-2 products.

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
# Sensor-specific L2 configuration
# ============================================================================

L2_SENSOR_CONFIG = {
    'OLI_TIRS': {
        'label': 'Landsat 8/9 OLI/TIRS',
        'sr_bands': [1, 2, 3, 4, 5, 6, 7],
        'st_band_suffix': 'ST_B10',
        'st_mtl_key': 'TEMPERATURE_MULT_BAND_ST_B10',
        'st_mtl_add_key': 'TEMPERATURE_ADD_BAND_ST_B10',
        'sr_wavelengths': [0.443, 0.482, 0.561, 0.655, 0.865, 1.609, 2.201],
        'st_wavelength': 10.895,
        'sr_names': [
            'Band 1 - Coastal Aerosol (0.443 um)', 'Band 2 - Blue (0.482 um)',
            'Band 3 - Green (0.561 um)', 'Band 4 - Red (0.655 um)',
            'Band 5 - NIR (0.865 um)', 'Band 6 - SWIR-1 (1.609 um)',
            'Band 7 - SWIR-2 (2.201 um)',
        ],
        'st_name': 'Surface Temperature (K) - Band 10 (10.895 um)',
    },
    'ETM': {
        'label': 'Landsat 7 ETM+',
        'sr_bands': [1, 2, 3, 4, 5, 7],
        'st_band_suffix': 'ST_B6',
        'st_mtl_key': 'TEMPERATURE_MULT_BAND_ST_B6',
        'st_mtl_add_key': 'TEMPERATURE_ADD_BAND_ST_B6',
        'sr_wavelengths': [0.485, 0.560, 0.660, 0.835, 1.650, 2.215],
        'st_wavelength': 11.450,
        'sr_names': [
            'Band 1 - Blue (0.485 um)', 'Band 2 - Green (0.560 um)',
            'Band 3 - Red (0.660 um)', 'Band 4 - NIR (0.835 um)',
            'Band 5 - SWIR-1 (1.650 um)', 'Band 7 - SWIR-2 (2.215 um)',
        ],
        'st_name': 'Surface Temperature (K) - Band 6 (11.450 um)',
    },
    'TM': {
        'label': 'Landsat 4/5 TM',
        'sr_bands': [1, 2, 3, 4, 5, 7],
        'st_band_suffix': 'ST_B6',
        'st_mtl_key': 'TEMPERATURE_MULT_BAND_ST_B6',
        'st_mtl_add_key': 'TEMPERATURE_ADD_BAND_ST_B6',
        'sr_wavelengths': [0.485, 0.560, 0.660, 0.830, 1.650, 2.215],
        'st_wavelength': 11.450,
        'sr_names': [
            'Band 1 - Blue (0.485 um)', 'Band 2 - Green (0.560 um)',
            'Band 3 - Red (0.660 um)', 'Band 4 - NIR (0.830 um)',
            'Band 5 - SWIR-1 (1.650 um)', 'Band 7 - SWIR-2 (2.215 um)',
        ],
        'st_name': 'Surface Temperature (K) - Band 6 (11.450 um)',
    },
}

# ST ancillary layer suffixes and names (same for all sensors)
ST_ANCILLARY_LAYERS = [
    ('ST_TRAD', 'Thermal Radiance'),
    ('ST_URAD', 'Upwelled Radiance'),
    ('ST_DRAD', 'Downwelled Radiance'),
    ('ST_ATRAN', 'Atmospheric Transmittance'),
    ('ST_EMIS', 'Emissivity'),
    ('ST_EMSD', 'Emissivity Stdev'),
    ('ST_CDIST', 'Cloud Distance'),
]

BAND_NAMES_ST_ANC = [name for _, name in ST_ANCILLARY_LAYERS]

ST_ANC_SCALE_FACTORS = {
    'ST_TRAD':  {'scale': 0.001, 'offset': 0.0},
    'ST_URAD':  {'scale': 0.001, 'offset': 0.0},
    'ST_DRAD':  {'scale': 0.001, 'offset': 0.0},
    'ST_ATRAN': {'scale': 0.0001, 'offset': 0.0},
    'ST_EMIS':  {'scale': 0.0001, 'offset': 0.0},
    'ST_EMSD':  {'scale': 0.0001, 'offset': 0.0},
    'ST_CDIST': {'scale': 0.01, 'offset': 0.0},
}


def detect_sensor(metadata):
    """Detect the Landsat sensor type from parsed MTL metadata."""
    root = metadata.get('LANDSAT_METADATA_FILE', metadata)
    img_attr = root.get('IMAGE_ATTRIBUTES', {})
    sensor_id = img_attr.get('SENSOR_ID', '')

    if sensor_id == 'TM':
        return 'TM'
    elif sensor_id in ('ETM', 'ETM+'):
        return 'ETM'
    elif sensor_id in ('OLI_TIRS', 'OLI', 'TIRS'):
        return 'OLI_TIRS'

    # Fallback from product ID
    prod_contents = root.get('PRODUCT_CONTENTS', {})
    product_id = prod_contents.get('LANDSAT_PRODUCT_ID', '')
    if product_id.startswith('LC0'):
        return 'OLI_TIRS'
    elif product_id.startswith('LE07'):
        return 'ETM'
    elif product_id.startswith('LT0'):
        return 'TM'

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


def extract_l2_calibration(metadata, config):
    """
    Extract L2SP calibration coefficients from parsed MTL metadata.
    Uses sensor config to find the right ST band key.
    """
    cal = {
        'sr_mult': {}, 'sr_add': {},
        'st_mult': None, 'st_add': None,
    }

    root = metadata.get('LANDSAT_METADATA_FILE', metadata)

    # Surface reflectance parameters
    sr_params = root.get('LEVEL2_SURFACE_REFLECTANCE_PARAMETERS', {})
    for key, value in sr_params.items():
        if not isinstance(value, (int, float)):
            continue
        m = re.match(r'REFLECTANCE_MULT_BAND_(\d+)', key)
        if m:
            cal['sr_mult'][int(m.group(1))] = value
        m = re.match(r'REFLECTANCE_ADD_BAND_(\d+)', key)
        if m:
            cal['sr_add'][int(m.group(1))] = value

    # Surface temperature parameters (sensor-specific key)
    st_params = root.get('LEVEL2_SURFACE_TEMPERATURE_PARAMETERS', {})
    cal['st_mult'] = st_params.get(config['st_mtl_key'])
    cal['st_add'] = st_params.get(config['st_mtl_add_key'])

    return cal


# ============================================================================
# Scene discovery
# ============================================================================

def discover_scenes(input_dir):
    """Scan for Landsat L2SP scenes (identified by _MTL.txt).
    When multiple MTL files exist in one directory, prefer L2SP/L2SR."""
    scenes = []

    def pick_l2_mtl(mtl_list):
        """From a list of MTL filenames, prefer L2SP/L2SR over L1."""
        l2 = [f for f in mtl_list if '_L2SP_' in f or '_L2SR_' in f]
        if l2:
            return l2[0]
        return mtl_list[0]

    mtl_files = [f for f in os.listdir(input_dir) if f.endswith('_MTL.txt')]
    if mtl_files:
        chosen = pick_l2_mtl(mtl_files)
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
                chosen = pick_l2_mtl(mtl_files)
                mtl_path = os.path.join(subdir, chosen)
                scene_id = chosen.replace('_MTL.txt', '')
                scenes.append((subdir, scene_id, mtl_path))

    return scenes


def find_sr_band_file(scene_dir, band_num):
    """Find surface reflectance band file (SR_B1 through SR_B7)."""
    suffix = f'_SR_B{band_num}.TIF'
    for f in os.listdir(scene_dir):
        if f.upper().endswith(suffix.upper()):
            return os.path.join(scene_dir, f)
    return None


def find_st_band_file(scene_dir, st_suffix):
    """Find surface temperature band file (ST_B10 or ST_B6)."""
    suffix = f'_{st_suffix}.TIF'
    for f in os.listdir(scene_dir):
        if f.upper().endswith(suffix.upper()):
            return os.path.join(scene_dir, f)
    return None


def find_st_ancillary_file(scene_dir, layer_suffix):
    """Find an ST ancillary file (ST_TRAD, ST_URAD, etc.)."""
    suffix = f'_{layer_suffix}.TIF'
    for f in os.listdir(scene_dir):
        if f.upper().endswith(suffix.upper()):
            return os.path.join(scene_dir, f)
    return None


# ============================================================================
# Band reading and calibration
# ============================================================================

def read_and_scale(filepath, mult, add, fill_value=0):
    """Read a band file and apply scale/offset: physical = MULT * DN + ADD."""
    with rasterio.open(filepath) as src:
        dn = src.read(1)
        profile = src.profile.copy()
    data = np.float32(dn)
    fill_mask = (dn == fill_value)
    data = mult * data + add
    data[fill_mask] = np.nan
    return data, profile


# ============================================================================
# ENVI header + output writers (same as L1 converter)
# ============================================================================

def build_envi_header(envi_path, height, width, num_bands, dtype, crs,
                      transform, wavelengths=None, band_names=None,
                      description='', interleave='bsq'):
    """Write a custom ENVI header file."""
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
            description=description, interleave='bsq')
        output_files.append(envi_dat)
        output_files.append(envi_hdr)

    return output_files


# ============================================================================
# Main processing
# ============================================================================

def process_scene(scene_dir, scene_id, mtl_path, output_dir,
                  write_geotiff=True, write_envi=True, overwrite=True):
    """Process a single Landsat L2SP scene. Auto-detects sensor."""
    all_outputs = []

    print(f"  Parsing MTL: {os.path.basename(mtl_path)}")
    metadata = parse_mtl(mtl_path)

    sensor_key = detect_sensor(metadata)
    config = L2_SENSOR_CONFIG[sensor_key]
    print(f"  Sensor: {config['label']}")

    cal = extract_l2_calibration(metadata, config)

    os.makedirs(output_dir, exist_ok=True)

    # --- Surface Reflectance ---
    print(f"\n  === Processing Surface Reflectance ===")
    sr_data = []
    ref_profile = None
    actual_names = []
    actual_wl = []

    for i, band_num in enumerate(config['sr_bands']):
        filepath = find_sr_band_file(scene_dir, band_num)
        if filepath is None:
            print(f"    WARNING: SR_B{band_num} not found, skipping.")
            continue

        mult = cal['sr_mult'].get(band_num)
        add = cal['sr_add'].get(band_num)
        if mult is None or add is None:
            print(f"    WARNING: No SR coefficients for Band {band_num}.")
            continue

        print(f"    Reading SR_B{band_num} (MULT={mult}, ADD={add})...")
        data, profile = read_and_scale(filepath, mult, add)
        sr_data.append(data)
        actual_names.append(config['sr_names'][i]
                            if i < len(config['sr_names'])
                            else f'SR Band {band_num}')
        actual_wl.append(config['sr_wavelengths'][i]
                         if i < len(config['sr_wavelengths'])
                         else 0.0)

        if ref_profile is None:
            ref_profile = profile
            print(f"    Image: {profile['width']} x {profile['height']}")

    if sr_data and ref_profile:
        sr_base = os.path.join(output_dir, f"{scene_id}_SR")
        if not overwrite and os.path.isfile(sr_base + '.tif'):
            print("    Output exists, skipping (overwrite=False)")
        else:
            outputs = write_multiband(
                sr_data, sr_base, ref_profile,
                band_names=actual_names, wavelengths=actual_wl,
                description=f'{config["label"]} L2SP Surface Reflectance '
                            f'({scene_id})',
                write_geotiff=write_geotiff, write_envi=write_envi)
            all_outputs.extend(outputs)

    # --- Surface Temperature ---
    print(f"\n  === Processing Surface Temperature ===")
    st_filepath = find_st_band_file(scene_dir, config['st_band_suffix'])

    if st_filepath and cal['st_mult'] is not None:
        print(f"    Reading {config['st_band_suffix']} "
              f"(MULT={cal['st_mult']}, ADD={cal['st_add']})...")
        st_data, st_profile = read_and_scale(
            st_filepath, cal['st_mult'], cal['st_add'])
        print(f"    Image: {st_profile['width']} x {st_profile['height']}")

        st_base = os.path.join(output_dir, f"{scene_id}_ST")
        if not overwrite and os.path.isfile(st_base + '.tif'):
            print("    Output exists, skipping (overwrite=False)")
        else:
            outputs = write_multiband(
                [st_data], st_base, st_profile,
                band_names=[config['st_name']],
                wavelengths=[config['st_wavelength']],
                description=f'{config["label"]} L2SP Surface Temperature '
                            f'in Kelvin ({scene_id})',
                write_geotiff=write_geotiff, write_envi=write_envi)
            all_outputs.extend(outputs)
    else:
        print(f"    WARNING: {config['st_band_suffix']} not found or "
              f"missing calibration.")

    # --- ST Ancillary ---
    print(f"\n  === Processing ST Ancillary Layers ===")
    anc_data = []
    anc_profile = None

    for layer_suffix, layer_name in ST_ANCILLARY_LAYERS:
        filepath = find_st_ancillary_file(scene_dir, layer_suffix)
        if filepath is None:
            print(f"    WARNING: {layer_suffix} not found, skipping.")
            continue

        sf = ST_ANC_SCALE_FACTORS.get(layer_suffix, {'scale': 1, 'offset': 0})
        print(f"    Reading {layer_suffix} "
              f"(scale={sf['scale']}, offset={sf['offset']})...")
        data, profile = read_and_scale(filepath, sf['scale'], sf['offset'])
        anc_data.append(data)
        if anc_profile is None:
            anc_profile = profile

    if anc_data and anc_profile:
        anc_base = os.path.join(output_dir, f"{scene_id}_ST_Ancillary")
        if not overwrite and os.path.isfile(anc_base + '.tif'):
            print("    Output exists, skipping (overwrite=False)")
        else:
            outputs = write_multiband(
                anc_data, anc_base, anc_profile,
                band_names=BAND_NAMES_ST_ANC[:len(anc_data)],
                wavelengths=None,
                description=f'{config["label"]} L2SP ST Ancillary ({scene_id})',
                write_geotiff=write_geotiff, write_envi=write_envi)
            all_outputs.extend(outputs)

    return all_outputs


def process_directory(input_dir, output_dir=None, write_geotiff=True,
                      write_envi=True, overwrite=True):
    """Process all Landsat L2SP scenes found in the input directory."""
    if output_dir is None:
        output_dir = os.path.join(input_dir, 'converted')

    print(f"Scanning for Landsat L2SP scenes in: {input_dir}")
    scenes = discover_scenes(input_dir)

    if not scenes:
        print("No Landsat L2SP scenes found.")
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
        print("Landsat Collection 2 Level-2 (L2SP) Converter")
        print("=" * 50)
        print(f"\nUsage: python {os.path.basename(__file__)} <input_dir> "
              f"[output_dir] [options]")
        print(f"\nSupported sensors (auto-detected from MTL):")
        print(f"  Landsat 8/9 OLI/TIRS  |  Landsat 7 ETM+  |  Landsat 4/5 TM")
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
