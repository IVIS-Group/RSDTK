"""
ECOSTRESS L1CT Radiance Converter
Merges individual L1CT GeoTIFF band files into multiband GeoTIFF and ENVI outputs.
Band order: radiance_1-5, cloud, water (7 bands total)

ECOSTRESS TIR Band Wavelengths:
  Band 1: 8.29 µm
  Band 2: 8.63 µm
  Band 3: 9.07 µm
  Band 4: 10.49 µm
  Band 5: 12.01 µm
"""

import os
import re
import glob
import numpy as np
import rasterio
from rasterio.transform import from_bounds
from collections import defaultdict

# ============================================================================
# ECOSTRESS L1CT filename pattern
# Example: ECOv002_L1CT_RAD_43484_016_50RQP_20260307T203953_0713_01_radiance_1.tif
# Components: sensor_version_level_orbit_scene_tile_datetime_build_iter_layername.tif
# ============================================================================

L1CT_PATTERN = re.compile(
    r'^(ECOv002_L1CT_RAD_\d+_\d+_\w+_\d{8}T\d{6}_\d{4}_\d{2})_(.+)\.tif$',
    re.IGNORECASE
)

# Band layers in desired output order
RADIANCE_BANDS = ['radiance_1', 'radiance_2', 'radiance_3', 'radiance_4', 'radiance_5']
MASK_BANDS = ['cloud', 'water']
ALL_BANDS = RADIANCE_BANDS + MASK_BANDS

# ECOSTRESS TIR wavelengths in micrometers (bands 1-5)
WAVELENGTHS_UM = [8.29, 8.63, 9.07, 10.49, 12.01]

# Band names for ENVI header
BAND_NAMES = [
    'Radiance Band 1 (8.29 um)',
    'Radiance Band 2 (8.63 um)',
    'Radiance Band 3 (9.07 um)',
    'Radiance Band 4 (10.49 um)',
    'Radiance Band 5 (12.01 um)',
    'Cloud Mask',
    'Water Mask'
]


def discover_granules(input_dir):
    """
    Scan input directory for L1CT GeoTIFF files and group them by granule base name.
    Returns dict: {granule_base: {layer_name: filepath, ...}, ...}
    """
    granules = defaultdict(dict)
    
    tif_files = glob.glob(os.path.join(input_dir, '*.tif'))
    
    for filepath in tif_files:
        filename = os.path.basename(filepath)
        match = L1CT_PATTERN.match(filename)
        if match:
            granule_base = match.group(1)
            layer_name = match.group(2).lower()
            granules[granule_base][layer_name] = filepath
    
    return dict(granules)


def validate_granule(granule_base, layer_files):
    """
    Check that all required bands are present for a granule.
    Returns (is_valid, missing_bands, available_bands)
    """
    available = set(layer_files.keys())
    required = set(ALL_BANDS)
    missing = required - available
    
    return len(missing) == 0, missing, available


def get_reference_profile(layer_files):
    """
    Read the rasterio profile from the first radiance band to use as reference
    for the output file.
    """
    ref_path = layer_files[RADIANCE_BANDS[0]]
    with rasterio.open(ref_path) as src:
        profile = src.profile.copy()
        height = src.height
        width = src.width
        crs = src.crs
        transform = src.transform
    
    return profile, height, width, crs, transform


def build_envi_header(envi_path, height, width, num_bands, dtype, crs, transform,
                      wavelengths=None, band_names=None, interleave='bsq'):
    """
    Write a custom ENVI header file with wavelength metadata.
    """
    # Map numpy/rasterio dtypes to ENVI data type codes
    dtype_map = {
        'uint8': 1,
        'int16': 2,
        'int32': 3,
        'float32': 4,
        'float64': 5,
        'uint16': 12,
        'uint32': 13,
        'int64': 14,
        'uint64': 15,
    }
    
    dtype_str = np.dtype(dtype).name
    envi_dtype = dtype_map.get(dtype_str, 4)  # Default to float32
    
    # Determine byte order (0 = little endian, 1 = big endian)
    byte_order = 0
    
    hdr_path = envi_path + '.hdr' if not envi_path.endswith('.hdr') else envi_path
    dat_path = hdr_path.replace('.hdr', '.dat')
    
    with open(hdr_path, 'w') as f:
        f.write('ENVI\n')
        f.write(f'description = {{ECOSTRESS L1CT Radiance Multiband Image}}\n')
        f.write(f'samples = {width}\n')
        f.write(f'lines = {height}\n')
        f.write(f'bands = {num_bands}\n')
        f.write(f'header offset = 0\n')
        f.write(f'file type = ENVI Standard\n')
        f.write(f'data type = {envi_dtype}\n')
        f.write(f'interleave = {interleave}\n')
        f.write(f'byte order = {byte_order}\n')
        
        # Map info: {proj_name, ref_x, ref_y, easting, northing, x_size, y_size, zone, datum}
        if transform is not None:
            x_size = abs(transform.a)
            y_size = abs(transform.e)
            ul_x = transform.c
            ul_y = transform.f
            f.write(f'map info = {{Geographic Lat/Lon, 1, 1, {ul_x}, {ul_y}, {x_size}, {y_size}, WGS-84}}\n')
        
        if crs is not None:
            f.write(f'coordinate system string = {{{crs.to_wkt()}}}\n')
        
        # Wavelengths for radiance bands (first 5 bands only)
        if wavelengths is not None:
            f.write('wavelength units = Micrometers\n')
            wl_str = ', '.join(f'{w:.2f}' for w in wavelengths)
            # Pad with 0.00 for mask bands
            n_mask = num_bands - len(wavelengths)
            if n_mask > 0:
                wl_str += ', ' + ', '.join(['0.00'] * n_mask)
            f.write(f'wavelength = {{\n  {wl_str}}}\n')
        
        if band_names is not None:
            bn_str = ',\n  '.join(band_names)
            f.write(f'band names = {{\n  {bn_str}}}\n')
    
    return hdr_path, dat_path


def convert_granule(granule_base, layer_files, output_dir, output_geotiff=True,
                    output_envi=True,
                    overwrite=True):
    """
    Convert a single L1CT granule to multiband GeoTIFF and/or ENVI format.
    
    Parameters:
        granule_base: str - Base granule name (e.g., ECOv002_L1CT_RAD_43484_016_...)
        layer_files: dict - {layer_name: filepath} mapping
        output_dir: str - Output directory path
        output_geotiff: bool - Whether to produce GeoTIFF output
        output_envi: bool - Whether to produce ENVI output
    
    Returns:
        list of output file paths
    """
    output_files = []
    
    # Validate
    is_valid, missing, available = validate_granule(granule_base, layer_files)
    if not is_valid:
        print(f"  WARNING: Missing bands for {granule_base}: {missing}")
        print(f"  Available: {available}")
        print(f"  Skipping this granule.")
        return output_files
    
    # Get reference spatial info from first radiance band
    profile, height, width, crs, transform = get_reference_profile(layer_files)
    
    print(f"  Image dimensions: {width} x {height}")
    print(f"  CRS: {crs}")
    
    # Read all bands into array
    num_bands = len(ALL_BANDS)
    
    # Determine output dtype - use float32 to accommodate both radiance and mask data
    out_dtype = np.float32
    multiband_data = np.zeros((num_bands, height, width), dtype=out_dtype)
    
    for i, band_name in enumerate(ALL_BANDS):
        filepath = layer_files[band_name]
        print(f"  Reading {band_name}...")
        with rasterio.open(filepath) as src:
            data = src.read(1)
            
            # Verify dimensions match
            if data.shape != (height, width):
                print(f"  WARNING: {band_name} shape {data.shape} doesn't match "
                      f"reference shape ({height}, {width}). Skipping granule.")
                return output_files
            
            multiband_data[i] = data.astype(out_dtype)
    
    # Create output directory if needed
    os.makedirs(output_dir, exist_ok=True)
    
    # --- GeoTIFF output ---
    if output_geotiff:
        gtiff_filename = f"{granule_base}_multiband.tif"
        gtiff_path = os.path.join(output_dir, gtiff_filename)
        
        gtiff_profile = {
            'driver': 'GTiff',
            'dtype': out_dtype,
            'width': width,
            'height': height,
            'count': num_bands,
            'crs': crs,
            'transform': transform,
            'nodata': np.nan,
            'compress': 'lzw',
            'tiled': True,
            'blockxsize': 256,
            'blockysize': 256,
        }
        
        print(f"  Writing GeoTIFF: {gtiff_filename}")
        with rasterio.open(gtiff_path, 'w', **gtiff_profile) as dst:
            for i in range(num_bands):
                dst.write(multiband_data[i], i + 1)
                dst.set_band_description(i + 1, BAND_NAMES[i])
        
        output_files.append(gtiff_path)
        print(f"  GeoTIFF written: {gtiff_path}")
    
    # --- ENVI output ---
    if output_envi:
        envi_basename = f"{granule_base}_multiband"
        envi_hdr_path = os.path.join(output_dir, envi_basename + '.hdr')
        envi_dat_path = os.path.join(output_dir, envi_basename + '.dat')
        
        print(f"  Writing ENVI: {envi_basename}.dat/.hdr")
        
        # Write binary data file (BSQ interleave: band sequential)
        multiband_data.tofile(envi_dat_path)
        
        # Write header
        build_envi_header(
            envi_hdr_path,
            height=height,
            width=width,
            num_bands=num_bands,
            dtype=out_dtype,
            crs=crs,
            transform=transform,
            wavelengths=WAVELENGTHS_UM,
            band_names=BAND_NAMES,
            interleave='bsq'
        )
        
        output_files.append(envi_dat_path)
        output_files.append(envi_hdr_path)
        print(f"  ENVI written: {envi_dat_path}")
    
    return output_files


def process_directory(input_dir, output_dir=None, output_geotiff=True, output_envi=True,
                      overwrite=True):
    """
    Process all L1CT granules found in the input directory.
    
    Parameters:
        input_dir: str - Directory containing L1CT GeoTIFF files
        output_dir: str - Output directory (defaults to input_dir/converted)
        output_geotiff: bool - Whether to produce GeoTIFF output
        output_envi: bool - Whether to produce ENVI output
    """
    if output_dir is None:
        output_dir = os.path.join(input_dir, 'converted')
    
    print(f"Scanning for ECOSTRESS L1CT files in: {input_dir}")
    granules = discover_granules(input_dir)
    
    if not granules:
        print("No ECOSTRESS L1CT granules found.")
        return
    
    print(f"Found {len(granules)} granule(s):\n")
    
    all_outputs = []
    for i, (granule_base, layer_files) in enumerate(sorted(granules.items()), 1):
        print(f"--- Granule {i}/{len(granules)}: {granule_base} ---")
        print(f"  Layers found: {sorted(layer_files.keys())}")
        
        outputs = convert_granule(
            granule_base, layer_files, output_dir,
            output_geotiff=output_geotiff,
            output_envi=output_envi,
            overwrite=overwrite
        )
        all_outputs.extend(outputs)
        print()
    
    print(f"Processing complete. {len(all_outputs)} output file(s) created in: {output_dir}")


# ============================================================================
# Main entry point
# ============================================================================

if __name__ == '__main__':
    import sys
    
    if len(sys.argv) < 2:
        print("ECOSTRESS L1CT Radiance Converter")
        print("=" * 40)
        print(f"\nUsage: python {os.path.basename(__file__)} <input_directory> [output_directory]")
        print(f"\nExample:")
        print(f"  python {os.path.basename(__file__)} D:\\Data\\ECOSTRESS\\L1CT")
        print(f"  python {os.path.basename(__file__)} D:\\Data\\ECOSTRESS\\L1CT D:\\Data\\ECOSTRESS\\output")
        print(f"\nThe script will scan the input directory for ECOSTRESS L1CT GeoTIFF files,")
        print(f"group them by granule, and merge them into multiband GeoTIFF and ENVI files.")
        print(f"\nOutput band order:")
        for i, name in enumerate(BAND_NAMES, 1):
            print(f"  Band {i}: {name}")
        sys.exit(0)
    
    input_dir = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else None
    
    if not os.path.isdir(input_dir):
        print(f"Error: Input directory does not exist: {input_dir}")
        sys.exit(1)
    
    process_directory(input_dir, output_dir)