# Decompiled with PyLingual (https://pylingual.io)
# Internal filename: 'AVIRIS_Converter.py'
# Bytecode version: 3.11a7e (3495)
# Source timestamp: 1970-01-01 00:00:00 UTC (0)

"""
AVIRIS-NG & Classic AVIRIS Converter
=====================================
Reads orthorectified AVIRIS-NG and Classic AVIRIS ENVI binary files
(L1B radiance and L2 reflectance) and outputs multiband GeoTIFF and/or
ENVI files with proper wavelength metadata, nodata handling, and optional
water vapor band exclusion.

Supported products:
  AVIRIS-NG L1B:  ang{YYYYMMDD}t{HHMMSS}_rdn_*_img      (radiance, µW/nm/cm²/sr)
  AVIRIS-NG L2:   ang{YYYYMMDD}t{HHMMSS}_corr_*_img      (surface reflectance)
  Classic L1B:    f{YYMMDD}t{NN}p{NN}r{NN}rdn_g_*ort_img (radiance)
  Classic L2:     f{YYMMDD}t{NN}p{NN}r{NN}_rfl            (reflectance)

All products are already orthorectified (UTM) with map info in the ENVI header.
The converter reads the ENVI binary+header, applies optional band filtering
(water vapor exclusion, bad bands list, wavelength subsetting), and writes
standardized GeoTIFF/ENVI output.

Dependencies: numpy, rasterio
"""
import os
import re
import sys
import math
import struct
import numpy as np
import warnings
try:
    import rasterio
    from rasterio.crs import CRS
    from rasterio.transform import Affine
except ImportError:
    print('ERROR: rasterio is required. Install with: pip install rasterio')
    sys.exit(1)
WATER_VAPOR_REGIONS = [(1334, 1431), (1788, 1980)]
ENVI_DTYPES = {1: (np.uint8, 'B', 1), 2: (np.int16, 'h', 2), 3: (np.int32, 'i', 4), 4: (np.float32, 'f', 4), 5: (np.float64, 'd', 8), 12: (np.uint16, 'H', 2), 13: (np.uint32, 'I', 4), 14: (np.int64, 'q', 8), 15: (np.uint64, 'Q', 8)}
def parse_envi_header(hdr_path):
    """
    Parse an ENVI .hdr file into a dictionary.
    Handles multi-line values enclosed in { }.
    """
    with open(hdr_path, 'r') as f:
        text = f.read()
    result = {}
    lines = text.split('\n')
    joined = []
    buf = ''
    in_braces = False
    for line in lines:
        if in_braces:
            buf += ' ' + line.strip()
            if '}' in line:
                in_braces = False
                joined.append(buf)
                buf = ''
        else:
            if '{' in line and '}' not in line:
                in_braces = True
                buf = line.strip()
            else:
                joined.append(line.strip())
    for line in joined:
        if '=' not in line:
            continue
        else:
            key, _, val = line.partition('=')
            key = key.strip().lower()
            val = val.strip()
            if val.startswith('{') and val.endswith('}'):
                    val = val[1:(-1)].strip()
            result[key] = val
    return result
def parse_float_list(s):
    """Parse a comma-separated string of floats."""
    if not s:
        return np.array([])
    else:
        return np.array([float(x.strip()) for x in s.split(',') if x.strip()])
def parse_map_info(map_info_str):
    """
    Parse ENVI map info string into components.
    Format: {proj, ref_x, ref_y, easting, northing, dx, dy, zone, hemi, datum, ...}
    May include rotation= parameter.
    """
    parts = [p.strip() for p in map_info_str.split(',')]
    info = {
        'projection': parts[0] if len(parts) > 0 else 'UTM',
        'ref_x': float(parts[1]) if len(parts) > 1 else 1.0,
        'ref_y': float(parts[2]) if len(parts) > 2 else 1.0,
        'easting': float(parts[3]) if len(parts) > 3 else 0.0,
        'northing': float(parts[4]) if len(parts) > 4 else 0.0,
        'dx': float(parts[5]) if len(parts) > 5 else 1.0,
        'dy': float(parts[6]) if len(parts) > 6 else 1.0,
        'zone': None,
        'hemisphere': 'North',
        'datum': 'WGS-84',
        'rotation': 0.0,
    }
    for i, part in enumerate(parts[7:], start=7):
        part_clean = part.strip()
        if part_clean.isdigit():
            info['zone'] = int(part_clean)
        else:
            if part_clean.lower() in ['north', 'south']:
                info['hemisphere'] = part_clean.capitalize()
            else:
                if 'wgs' in part_clean.lower():
                    info['datum'] = part_clean
                else:
                    if 'rotation' in part_clean.lower():
                        rot_match = re.search('rotation\\s*=\\s*([-\\d.]+)', part_clean)
                        if rot_match:
                            info['rotation'] = float(rot_match.group(1))
    return info
def build_transform(map_info):
    """Build a rasterio Affine transform from parsed map info."""
    ulx = map_info['easting'] - (map_info['ref_x'] - 1) * map_info['dx']
    uly = map_info['northing'] + (map_info['ref_y'] - 1) * map_info['dy']
    dx = map_info['dx']
    dy = map_info['dy']
    rotation = map_info.get('rotation', 0.0)
    if abs(rotation) > 1e-06:
        rad = math.radians(rotation)
        cos_r = math.cos(rad)
        sin_r = math.sin(rad)
        transform = Affine(dx * cos_r, -dx * sin_r, ulx, -dy * sin_r, -dy * cos_r, uly)
    else:
        transform = Affine(dx, 0.0, ulx, 0.0, -dy, uly)
    return transform
def build_crs(map_info):
    """Build a rasterio CRS from parsed map info using WKT to avoid PROJ db issues."""
    zone = map_info.get('zone')
    hemi = map_info.get('hemisphere', 'North')
    if zone and map_info['projection'].upper() == 'UTM':
        south = hemi.lower() == 'south'
        false_northing = 10000000.0 if south else 0.0
        central_meridian = (zone - 1) * 6 - 180 + 3
        hemi_label = 'S' if south else 'N'
        epsg = 32700 + zone if south else 32600 + zone
        wkt = f'PROJCS[\"WGS 84 / UTM zone {zone}{hemi_label}\",GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",{central_meridian}],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",{false_northing}],UNIT[\"metre\",1],AUTHORITY[\"EPSG\",\"{epsg}\"]]'
        return CRS.from_wkt(wkt)
    else:
        wkt = 'GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433],AUTHORITY[\"EPSG\",\"4326\"]]'
        return CRS.from_wkt(wkt)
def read_envi_band(filepath, hdr, band_idx):
    """
    Read a single band from an ENVI binary file.

    Parameters:
        filepath: path to binary data file
        hdr: parsed header dict
        band_idx: 0-based band index

    Returns:
        2D numpy array (lines, samples)
    """
    samples = int(hdr['samples'])
    lines = int(hdr['lines'])
    bands = int(hdr['bands'])
    dtype_code = int(hdr['data type'])
    interleave = hdr.get('interleave', 'bil').lower()
    offset = int(hdr.get('header offset', 0))
    byte_order = int(hdr.get('byte order', 0))
    np_dtype = ENVI_DTYPES[dtype_code][0]
    if byte_order == 0:
        dt = np.dtype(np_dtype).newbyteorder('<')
    else:
        dt = np.dtype(np_dtype).newbyteorder('>')
    pixel_size = dt.itemsize
    if interleave == 'bsq':
        band_offset = offset + band_idx * lines * samples * pixel_size
        data = np.fromfile(filepath, dtype=dt, count=lines * samples, offset=band_offset)
        return data.reshape(lines, samples)
    else:
        if interleave == 'bil':
            data = np.empty((lines, samples), dtype=dt)
            line_size = bands * samples * pixel_size
            for row in range(lines):
                row_offset = offset + row * line_size + band_idx * samples * pixel_size
                data[row, :] = np.fromfile(filepath, dtype=dt, count=samples, offset=row_offset)
            return data
        else:
            if interleave == 'bip':
                data = np.empty((lines, samples), dtype=dt)
                line_size = samples * bands * pixel_size
                for row in range(lines):
                    row_data = np.fromfile(filepath, dtype=dt, count=samples * bands, offset=offset + row * line_size)
                    data[row, :] = row_data.reshape(samples, bands)[:, band_idx]
                return data
            else:
                raise ValueError(f'Unknown interleave: {interleave}')
def read_envi_all_bands(filepath, hdr, band_indices=None):
    """
    Read multiple bands from an ENVI binary file using memory mapping.

    Parameters:
        filepath: path to binary data file
        hdr: parsed header dict
        band_indices: list of 0-based band indices, or None for all

    Returns:
        3D numpy array (bands, lines, samples)
    """
    samples = int(hdr['samples'])
    lines = int(hdr['lines'])
    bands = int(hdr['bands'])
    dtype_code = int(hdr['data type'])
    interleave = hdr.get('interleave', 'bil').lower()
    offset = int(hdr.get('header offset', 0))
    byte_order = int(hdr.get('byte order', 0))
    np_dtype = ENVI_DTYPES[dtype_code][0]
    if byte_order == 0:
        dt = np.dtype(np_dtype).newbyteorder('<')
    else:
        dt = np.dtype(np_dtype).newbyteorder('>')
    if band_indices is None:
        band_indices = list(range(bands))
    n_out = len(band_indices)
    if interleave == 'bsq':
        mmap = np.memmap(filepath, dtype=dt, mode='r', offset=offset, shape=(bands, lines, samples))
        result = np.empty((n_out, lines, samples), dtype=np.float32)
        for i, bi in enumerate(band_indices):
            result[i] = mmap[bi].astype(np.float32)
    else:
        if interleave == 'bil':
            mmap = np.memmap(filepath, dtype=dt, mode='r', offset=offset, shape=(lines, bands, samples))
            result = np.empty((n_out, lines, samples), dtype=np.float32)
            for i, bi in enumerate(band_indices):
                result[i] = mmap[:, bi, :].astype(np.float32)
        else:
            if interleave == 'bip':
                mmap = np.memmap(filepath, dtype=dt, mode='r', offset=offset, shape=(lines, samples, bands))
                result = np.empty((n_out, lines, samples), dtype=np.float32)
                for i, bi in enumerate(band_indices):
                    result[i] = mmap[:, :, bi].astype(np.float32)
            else:
                raise ValueError(f'Unknown interleave: {interleave}')
    del mmap
    return result
NG_RDN_PATTERN = re.compile('^(ang\\d{8}t\\d{6})_rdn_([^.]+)_img$', re.IGNORECASE)
NG_CORR_PATTERN = re.compile('^(ang\\d{8}t\\d{6})_corr_([^.]+)_img$', re.IGNORECASE)
CLASSIC_RDN_ORT_PATTERN = re.compile('^(f\\d{6}t\\d{2}p\\d{2}r\\d{2})rdn_g(?:_sc\\d+)?_ort_img$', re.IGNORECASE)
CLASSIC_RFL_PATTERN = re.compile('^(f\\d{6}t\\d{2}p\\d{2}r\\d{2})_rfl$', re.IGNORECASE)
def discover_files(input_dir):
    """
    Scan for AVIRIS-NG and Classic AVIRIS data files.

    Returns list of tuples: (filepath, hdr_path, flight_id, product_type, variant)
    where variant is \'ng\' or \'classic\' and product_type is \'rdn\' or \'rfl\'/\'corr\'.
    """
    files = []
    def scan_dir(d):
        try:
            entries = os.listdir(d)
        except OSError:
            return None
        for fname in sorted(entries):
            fpath = os.path.join(d, fname)
            if not os.path.isfile(fpath):
                continue
            else:
                if fname.lower().endswith(('.hdr', '.glt', '.igm', '.loc', '.obs', '.obs_ort', '.gain', '.rccf2', '.readme', '.eph', '.lonlat_eph', '.plog')):
                    continue
                else:
                    m = NG_RDN_PATTERN.match(fname)
                    if m:
                        hdr = fpath + '.hdr'
                        if os.path.isfile(hdr):
                            files.append((fpath, hdr, m.group(1), 'rdn', 'ng'))
                        continue
                    else:
                        m = NG_CORR_PATTERN.match(fname)
                        if m:
                            hdr = fpath + '.hdr'
                            if os.path.isfile(hdr):
                                files.append((fpath, hdr, m.group(1), 'corr', 'ng'))
                            continue
                        else:
                            m = CLASSIC_RDN_ORT_PATTERN.match(fname)
                            if m:
                                hdr = fpath + '.hdr'
                                if os.path.isfile(hdr):
                                    files.append((fpath, hdr, m.group(1), 'rdn', 'classic'))
                                continue
                            else:
                                m = CLASSIC_RFL_PATTERN.match(fname)
                                if m:
                                    hdr = fpath + '.hdr'
                                    if os.path.isfile(hdr):
                                        files.append((fpath, hdr, m.group(1), 'rfl', 'classic'))
                                    continue
    scan_dir(input_dir)
    for entry in os.listdir(input_dir):
        subdir = os.path.join(input_dir, entry)
        if os.path.isdir(subdir):
            scan_dir(subdir)
    return sorted(files, key=lambda x: (x[4], x[2], x[3]))
def get_band_mask(wavelengths, bbl=None, exclude_water=True, wl_min=None, wl_max=None):
    """
    Build a boolean mask of bands to include.

    Parameters:
        wavelengths: array of band center wavelengths (nm)
        bbl: bad bands list (1=good, 0=bad), or None
        exclude_water: auto-exclude water vapor absorption bands
        wl_min: minimum wavelength to include (nm), or None
        wl_max: maximum wavelength to include (nm), or None

    Returns:
        boolean array, True = include
    """
    n = len(wavelengths)
    mask = np.ones(n, dtype=bool)
    if bbl is not None and len(bbl) == n:
            mask &= bbl > 0.5
    if exclude_water:
        for wl_lo, wl_hi in WATER_VAPOR_REGIONS:
            mask &= ~((wavelengths >= wl_lo) & (wavelengths <= wl_hi))
    if wl_min is not None:
        mask &= wavelengths >= wl_min
    if wl_max is not None:
        mask &= wavelengths <= wl_max
    return mask
def write_geotiff(filepath, data, transform, crs, wavelengths=None, nodata=np.nan):
    """
    Write a multiband GeoTIFF with DEFLATE compression.

    Parameters:
        filepath: output path
        data: 3D array (bands, lines, samples), float32
        transform: rasterio Affine
        crs: rasterio CRS
        wavelengths: array of wavelengths for band descriptions
        nodata: nodata value
    """
    n_bands, height, width = data.shape
    profile = {'driver': 'GTiff', 'dtype': 'float32', 'width': width, 'height': height, 'count': n_bands, 'crs': crs, 'transform': transform, 'nodata': nodata, 'compress': 'deflate', 'tiled': True, 'blockxsize': 256, 'blockysize': 256, 'predictor': 2}
    with rasterio.open(filepath, 'w', **profile) as dst:
        for i in range(n_bands):
            dst.write(data[i], i + 1)
            if wavelengths is not None and i < len(wavelengths):
                    dst.set_band_description(i + 1, f'Band {i + 1} ({wavelengths[i]:.2f} nm)')
def write_envi(filepath, data, transform, crs, wavelengths, fwhm=None, nodata_val=(-9999)):
    """
    Write ENVI binary (BSQ) + header with wavelength metadata.
    Uses manual header writing for reliable wavelength embedding.

    Parameters:
        filepath: output path (no extension — writes .dat and .hdr)
        data: 3D array (bands, lines, samples), float32
        transform: rasterio Affine
        crs: rasterio CRS
        wavelengths: array of wavelengths (nm)
        fwhm: array of FWHM values (nm), or None
        nodata_val: nodata value for ENVI output
    """
    dat_path = filepath + '.dat'
    hdr_path = filepath + '.hdr'
    n_bands, height, width = data.shape
    out_data = data.copy()
    out_data[np.isnan(out_data)] = nodata_val
    out_data.astype('<f4').tofile(dat_path)
    map_info_str = _build_map_info_string(transform, crs)
    crs_wkt = crs.to_wkt() if crs else ''
    hdr_lines = ['ENVI', f'description = {{{dat_path}}}', f'samples = {width}', f'lines   = {height}', f'bands   = {n_bands}', 'header offset = 0', 'file type = ENVI Standard', 'data type = 4', 'interleave = bsq', 'byte order = 0']
    if map_info_str:
        hdr_lines.append(f'map info = {{{map_info_str}}}')
    if crs_wkt:
        hdr_lines.append(f'coordinate system string = {{{crs_wkt}}}')
    if wavelengths is not None and len(wavelengths) > 0:
            wl_str = ', '.join((f'{w:.4f}' for w in wavelengths))
            hdr_lines.append(f'wavelength = {{{wl_str}}}')
            hdr_lines.append('wavelength units = Nanometers')
    if fwhm is not None and len(fwhm) > 0:
            fwhm_str = ', '.join((f'{f:.4f}' for f in fwhm))
            hdr_lines.append(f'fwhm = {{{fwhm_str}}}')
    band_names = ', '.join((f'Band {i + 1} ({wavelengths[i]:.2f} nm)' for i in range(n_bands)))
    hdr_lines.append(f'band names = {{{band_names}}}')
    hdr_lines.append(f'data ignore value = {nodata_val}')
    with open(hdr_path, 'w') as f:
        f.write('\n'.join(hdr_lines) + '\n')
def _build_map_info_string(transform, crs):
    """Build an ENVI map info string from transform and CRS."""
    if crs is None:
        return ''
    crs_dict = crs.to_dict()
    epsg = None
    try:
        epsg = crs.to_epsg()
    except:
        pass
    if epsg and 32601 <= epsg <= 32660:
        zone = epsg - 32600
        hemi = 'North'
    elif epsg and 32701 <= epsg <= 32760:
        zone = epsg - 32700
        hemi = 'South'
    else:
        zone = None
        hemi = 'North'
    ulx = transform.c
    uly = transform.f
    dx = abs(transform.a)
    dy = abs(transform.e)
    if zone:
        return f'UTM, 1, 1, {ulx}, {uly}, {dx}, {dy}, {zone}, {hemi}, WGS-84, units=Meters'
    else:
        return f'Geographic Lat/Lon, 1, 1, {ulx}, {uly}, {dx}, {dy}, WGS-84'
def process_file(filepath, hdr_path, flight_id, product_type, variant, output_dir, geotiff=True, envi=False, exclude_water=True, use_bbl=True, wl_min=None, wl_max=None, overwrite=True, target_resolution=None):
    """
    Process a single AVIRIS file.

    Parameters:
        filepath: path to ENVI binary file
        hdr_path: path to companion .hdr file
        flight_id: extracted flight line identifier
        product_type: \'rdn\', \'corr\', or \'rfl\'
        variant: \'ng\' or \'classic\'
        output_dir: root output directory
        geotiff: write GeoTIFF output
        envi: write ENVI output
        exclude_water: auto-exclude water vapor bands
        use_bbl: use bad bands list from header (NG only)
        wl_min: minimum wavelength filter (nm)
        wl_max: maximum wavelength filter (nm)
        overwrite: overwrite existing files
        target_resolution: optional target pixel size in meters for resampling
    """
    variant_label = 'AVIRIS-NG' if variant == 'ng' else 'Classic AVIRIS'
    prod_label = {'rdn': 'L1B Radiance', 'corr': 'L2 Reflectance', 'rfl': 'L2 Reflectance'}[product_type]
    print(f"\n{'======================================================================'}")
    print(f'Processing {variant_label} {prod_label}: {flight_id}')
    print(f'  File: {os.path.basename(filepath)}')
    hdr = parse_envi_header(hdr_path)
    samples = int(hdr['samples'])
    lines = int(hdr['lines'])
    bands = int(hdr['bands'])
    print(f'  Dimensions: {lines} lines × {samples} samples × {bands} bands')
    wavelengths = parse_float_list(hdr.get('wavelength', ''))
    if len(wavelengths) != bands:
        print(f'  WARNING: wavelength count ({len(wavelengths)}) != band count ({bands}). Skipping file.')
        return
    else:
        wl_units = hdr.get('wavelength units', 'Nanometers').lower()
        if 'micro' in wl_units:
            wavelengths *= 1000.0
            print('  Converted wavelengths from µm to nm')
        print(f'  Wavelength range: {wavelengths.min():.1f} - {wavelengths.max():.1f} nm')
        fwhm = parse_float_list(hdr.get('fwhm', ''))
        if len(fwhm) != bands:
            fwhm = None
        else:
            if 'micro' in wl_units and fwhm is not None:
                    fwhm *= 1000.0
        bbl = None
        if use_bbl and 'bbl' in hdr:
                bbl = parse_float_list(hdr['bbl'])
                if len(bbl) != bands:
                    bbl = None
                else:
                    n_bad = int(np.sum(bbl < 0.5))
                    if n_bad > 0:
                        print(f'  Bad bands list: {n_bad} bad bands flagged')
        band_mask = get_band_mask(wavelengths, bbl=bbl, exclude_water=exclude_water, wl_min=wl_min, wl_max=wl_max)
        band_indices = np.where(band_mask)[0]
        n_out = len(band_indices)
        if n_out == 0:
            print('  WARNING: No bands remaining after filtering. Skipping.')
            return
        else:
            excluded = bands - n_out
            if excluded > 0:
                print(f'  Band filtering: {n_out} of {bands} bands selected ({excluded} excluded)')
            map_info_str = hdr.get('map info', '')
            if not map_info_str:
                print('  WARNING: No map info found in header. Skipping file.')
                return
            else:
                map_info = parse_map_info(map_info_str)
                transform = build_transform(map_info)
                crs = build_crs(map_info)
                print(f"  CRS: EPSG:{(crs.to_epsg() if crs.to_epsg() else 'unknown')} | Pixel: {map_info['dx']:.1f} m | Rotation: {map_info['rotation']:.2f}°")
                prod_suffix = {'rdn': 'RDN', 'corr': 'RFL', 'rfl': 'RFL'}[product_type]
                prefix = 'AVIRIS-NG' if variant == 'ng' else 'AVIRIS'
                out_base = f'{prefix}_{flight_id}_{prod_suffix}'
                tiff_dir = os.path.join(output_dir, 'GeoTIFF')
                envi_dir = os.path.join(output_dir, 'ENVI')
                tiff_path = os.path.join(tiff_dir, out_base + '.tif')
                envi_path = os.path.join(envi_dir, out_base)
                if not overwrite:
                    tiff_exists = os.path.isfile(tiff_path) if geotiff else False
                    envi_exists = os.path.isfile(envi_path + '.dat') if envi else False
                    if geotiff and tiff_exists and (not envi or envi_exists):
                        print('  Output exists, skipping (overwrite=False)')
                        return
                    else:
                        if envi and envi_exists and (not geotiff or tiff_exists):
                            print('  Output exists, skipping (overwrite=False)')
                            return
                print(f'  Reading {n_out} bands...')
                data = read_envi_all_bands(filepath, hdr, band_indices=band_indices.tolist())
                out_wavelengths = wavelengths[band_indices]
                out_fwhm = fwhm[band_indices] if fwhm is not None else None
                nodata_str = hdr.get('data ignore value', '')
                if nodata_str:
                    try:
                        src_nodata = float(nodata_str)
                    except ValueError:
                        src_nodata = (-9999.0)
                else:
                    src_nodata = (-9999.0)
                if not np.isnan(src_nodata):
                    nodata_mask = np.isclose(data, src_nodata)
                    data[nodata_mask] = np.nan
                if target_resolution is not None and target_resolution > map_info['dx']:
                        data, transform = _resample_block_average(data, transform, map_info['dx'], target_resolution)
                        lines, samples = (data.shape[1], data.shape[2])
                        print(f'  Resampled to {target_resolution:.1f} m ({lines} × {samples} pixels)')
                if geotiff:
                    os.makedirs(tiff_dir, exist_ok=True)
                    print(f'  Writing GeoTIFF: {out_base}.tif ({n_out} bands)')
                    write_geotiff(tiff_path, data, transform, crs, wavelengths=out_wavelengths, nodata=np.nan)
                if envi:
                    os.makedirs(envi_dir, exist_ok=True)
                    print(f'  Writing ENVI: {out_base}.dat/.hdr ({n_out} bands)')
                    write_envi(envi_path, data, transform, crs, wavelengths=out_wavelengths, fwhm=out_fwhm, nodata_val=(-9999))
                print(f'  Done: {flight_id}')
def _resample_block_average(data, transform, src_res, target_res):
    """
    Resample by block averaging to a coarser resolution.

    Parameters:
        data: 3D array (bands, lines, samples)
        transform: source Affine transform
        src_res: source pixel size (m)
        target_res: target pixel size (m)

    Returns:
        (resampled_data, new_transform)
    """
    factor = int(round(target_res / src_res))
    if factor <= 1:
        return (data, transform)
    else:
        n_bands, lines, samples = data.shape
        new_lines = lines // factor
        new_samples = samples // factor
        trimmed = data[:, :new_lines * factor, :new_samples * factor]
        reshaped = trimmed.reshape(n_bands, new_lines, factor, new_samples, factor)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            result = np.nanmean(reshaped, axis=(2, 4)).astype(np.float32)
        new_transform = Affine(transform.a * factor, transform.b * factor, transform.c, transform.d * factor, transform.e * factor, transform.f)
        return (result, new_transform)
def process_directory(input_dir, output_dir=None, geotiff=True, envi=False, exclude_water=True, use_bbl=True, wl_min=None, wl_max=None, overwrite=True, target_resolution=None, bbox=None, spectral_opts=None):
    """
    Batch process AVIRIS-NG and Classic AVIRIS files.

    Parameters:
        input_dir: directory containing AVIRIS data
        output_dir: output directory (default: input_dir/AVIRIS_output)
        geotiff: write GeoTIFF output
        envi: write ENVI output
        exclude_water: auto-exclude water vapor absorption bands
        use_bbl: use bad bands list from NG headers
        wl_min: minimum wavelength filter (nm)
        wl_max: maximum wavelength filter (nm)
        overwrite: overwrite existing output files
        target_resolution: optional target pixel size for resampling (m)
        bbox: bounding box tuple (min_lon, min_lat, max_lon, max_lat) — reserved
        spectral_opts: spectral filtering options dict — reserved
    """
    if output_dir is None:
        output_dir = os.path.join(input_dir, 'AVIRIS_output')
    if not geotiff and (not envi):
            geotiff = True
    print('AVIRIS-NG / Classic AVIRIS Converter')
    print(f"{'======================================================================'}")
    print(f'Input:  {input_dir}')
    print(f'Output: {output_dir}')
    print(f"Format: {('GeoTIFF' if geotiff else '')}{('/' if geotiff and envi else '')}{('ENVI' if envi else '')}")
    options_parts = []
    if exclude_water:
        options_parts.append('water vapor exclusion')
    if use_bbl:
        options_parts.append('BBL filtering (NG)')
    if wl_min or wl_max:
        wl_range = f"{wl_min or '...'}-{wl_max or '...'} nm"
        options_parts.append(f'wavelength range: {wl_range}')
    if target_resolution:
        options_parts.append(f'resample to {target_resolution}m')
    if options_parts:
        print(f"Options: {', '.join(options_parts)}")
    files = discover_files(input_dir)
    if not files:
        print('\nNo AVIRIS-NG or Classic AVIRIS files found.')
        print('Expected filename patterns:')
        print('  NG L1B:      ang{YYYYMMDD}t{HHMMSS}_rdn_*_img + .hdr')
        print('  NG L2:       ang{YYYYMMDD}t{HHMMSS}_corr_*_img + .hdr')
        print('  Classic L1B: f{YYMMDD}t{NN}p{NN}r{NN}rdn_g_*ort_img + .hdr')
        print('  Classic L2:  f{YYMMDD}t{NN}p{NN}r{NN}_rfl + .hdr')
        return
    else:
        ng_count = sum((1 for f in files if f[4] == 'ng'))
        classic_count = sum((1 for f in files if f[4] == 'classic'))
        print(f'\nFound {len(files)} file(s): {ng_count} AVIRIS-NG, {classic_count} Classic AVIRIS')
        for i, (fpath, hdr_path, fid, ptype, var) in enumerate(files, 1):
            print(f'\n[{i}/{len(files)}]', end='')
            try:
                process_file(fpath, hdr_path, fid, ptype, var, output_dir=output_dir, geotiff=geotiff, envi=envi, exclude_water=exclude_water, use_bbl=use_bbl, wl_min=wl_min, wl_max=wl_max, overwrite=overwrite, target_resolution=target_resolution)
            except Exception as e:
                print(f'\n  ERROR processing {fid}: {e}')
                import traceback
                traceback.print_exc()
        print(f"\n{'======================================================================'}")
        print(f'Processing complete. Output: {output_dir}')
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='AVIRIS-NG & Classic AVIRIS Converter — Convert orthorectified ENVI data to GeoTIFF/ENVI')
    parser.add_argument('input_dir', help='Input directory with AVIRIS data')
    parser.add_argument('-o', '--output', help='Output directory')
    parser.add_argument('--geotiff', action='store_true', default=True, help='Write GeoTIFF (default)')
    parser.add_argument('--envi', action='store_true', help='Write ENVI output')
    parser.add_argument('--both', action='store_true', help='Write both GeoTIFF and ENVI')
    parser.add_argument('--no-water-exclude', action='store_true', help='Keep water vapor absorption bands')
    parser.add_argument('--no-bbl', action='store_true', help='Ignore bad bands list')
    parser.add_argument('--wl-min', type=float, help='Minimum wavelength (nm)')
    parser.add_argument('--wl-max', type=float, help='Maximum wavelength (nm)')
    parser.add_argument('--resolution', type=float, help='Target resolution for resampling (m)')
    parser.add_argument('--no-overwrite', action='store_true', help='Skip existing output files')
    args = parser.parse_args()
    do_geotiff = True
    do_envi = args.envi
    if args.both:
        do_geotiff = True
        do_envi = True
    process_directory(args.input_dir, output_dir=args.output, geotiff=do_geotiff, envi=do_envi, exclude_water=not args.no_water_exclude, use_bbl=not args.no_bbl, wl_min=args.wl_min, wl_max=args.wl_max, overwrite=not args.no_overwrite, target_resolution=args.resolution)