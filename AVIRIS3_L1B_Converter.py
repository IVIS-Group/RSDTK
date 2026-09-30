"""
AVIRIS-3/5 L1B Calibrated Radiance Converter
==============================================
Reads AVIRIS-3 and AVIRIS-5 L1B NetCDF4 files, applies GLT-based
orthorectification (UTM projection), and outputs multiband GeoTIFF
and/or ENVI files.

AVIRIS-3: ~284 bands, 390-2500 nm (~7.4 nm spacing), VSWIR
AVIRIS-5: ~424 bands, 390-2500 nm (~5 nm spacing), VSWIR
Spatial resolution: altitude-dependent (typically 4-20m)
Uses GLT (Geometric Lookup Table) for fast orthorectification.

Dependencies: h5py, numpy, rasterio
"""

import os
import re
import sys
import numpy as np
import h5py
import warnings

try:
    import rasterio
    from rasterio.crs import CRS
    from rasterio.transform import Affine
except ImportError:
    print("ERROR: rasterio is required. Install with: pip install rasterio")
    sys.exit(1)


# ============================================================================
# Constants
# ============================================================================

# Filename patterns — accepts both AV3 and AV5
# AV320230927t190252_001_L1B_RDN_V01.nc  (AVIRIS-3, facility collection)
# AV320230716t204543_L1B_RDN_001_V01.nc  (AVIRIS-3, ABoVE collection)
# AV520250603t211440_012_L1B_RDN_5b9e1cc2_RDN.nc  (AVIRIS-5)
AV_L1B_PATTERN = re.compile(
    r'^(AV[35]\d{8}t\d{6}).*L1B.*RDN.*\.nc$',
    re.IGNORECASE
)

# Water vapor absorption regions to auto-exclude (nm)
WATER_VAPOR_REGIONS = [
    (1334, 1440),   # 1.4 µm H2O (extended to catch edge channels)
    (1788, 1980),   # 1.9 µm H2O
]


# ============================================================================
# File discovery
# ============================================================================

def discover_files(input_dir):
    """Scan for AVIRIS-3/5 L1B NetCDF files."""
    files = []

    def scan_dir(d):
        try:
            entries = os.listdir(d)
        except OSError:
            return
        for f in entries:
            if AV_L1B_PATTERN.match(f):
                match = re.match(r'^(AV[35]\d{8}t\d{6})', f)
                gid = match.group(1) if match else f.replace('.nc', '')
                files.append((os.path.join(d, f), gid))

    scan_dir(input_dir)
    for entry in os.listdir(input_dir):
        subdir = os.path.join(input_dir, entry)
        if os.path.isdir(subdir):
            scan_dir(subdir)

    return sorted(files, key=lambda x: x[1])


# ============================================================================
# Data reading
# ============================================================================

def read_av_data(filepath, exclude_water=True, spectral_opts=None, bbox=None):
    """
    Read AVIRIS-3/5 L1B radiance, wavelengths, GLT, and CRS info.

    Parameters:
        filepath: path to NetCDF file
        exclude_water: auto-exclude water vapor absorption bands
        spectral_opts: dict with spectral filtering options
        bbox: bounding box tuple (min_lon, min_lat, max_lon, max_lat)

    Returns dict with keys:
        radiance, wavelengths, fwhm, glt_line, glt_sample,
        transform, crs, sensor_name
    """
    print(f"  Reading: {os.path.basename(filepath)}")

    with h5py.File(filepath, 'r') as f:
        # Determine sensor from filename or attributes
        fname = os.path.basename(filepath)
        if fname.upper().startswith('AV5'):
            sensor_name = 'AVIRIS-5'
        else:
            sensor_name = 'AVIRIS-3'

        # Read radiance data — bands-first: (bands, lines, samples)
        rad = f['radiance/radiance'][:]
        wavelengths = f['radiance/wavelength'][:]
        fwhm_arr = f['radiance/fwhm'][:]

        n_bands, n_lines, n_samples = rad.shape
        print(f"  {sensor_name}: {n_bands} bands, {n_lines} lines × "
              f"{n_samples} samples")
        print(f"  Wavelengths: {wavelengths.min():.1f} - "
              f"{wavelengths.max():.1f} nm")

        # Read GLT
        glt_line = f['geolocation_lookup_table/line'][:]
        glt_sample = f['geolocation_lookup_table/sample'][:]
        easting = f['geolocation_lookup_table/easting'][:]
        northing = f['geolocation_lookup_table/northing'][:]

        # Read CRS from transverse_mercator attributes
        tm_attrs = dict(f['transverse_mercator'].attrs)

        # Extract geotransform
        if 'GeoTransform' in tm_attrs:
            gt_str = tm_attrs['GeoTransform']
            if isinstance(gt_str, bytes):
                gt_str = gt_str.decode('utf-8')
            gt = [float(x) for x in gt_str.split()]
            transform = Affine(gt[1], gt[2], gt[0], gt[4], gt[5], gt[3])
        else:
            # Build from easting/northing arrays
            dx = float(np.diff(easting[:2])[0]) if len(easting) > 1 else 1.0
            dy = float(np.diff(northing[:2])[0]) if len(northing) > 1 else -1.0
            transform = Affine(dx, 0.0, float(easting[0]),
                               0.0, dy, float(northing[0]))

        # Build CRS
        if 'spatial_ref' in tm_attrs:
            wkt = tm_attrs['spatial_ref']
            if isinstance(wkt, bytes):
                wkt = wkt.decode('utf-8')
            crs = CRS.from_wkt(wkt)
        elif 'crs_wkt' in tm_attrs:
            wkt = tm_attrs['crs_wkt']
            if isinstance(wkt, bytes):
                wkt = wkt.decode('utf-8')
            crs = CRS.from_wkt(wkt)
        else:
            # Fallback: try to determine UTM zone from easting/northing
            zone = int(tm_attrs.get('utm_zone_number',
                       tm_attrs.get('zone', 0)))
            if zone > 0:
                crs = CRS.from_epsg(32600 + zone)
            else:
                print("  WARNING: Could not determine CRS, defaulting to "
                      "EPSG:32610")
                crs = CRS.from_epsg(32610)

    # --- Band filtering ---
    band_mask = np.ones(n_bands, dtype=bool)

    # Water vapor exclusion
    if exclude_water:
        for wl_lo, wl_hi in WATER_VAPOR_REGIONS:
            in_region = (wavelengths >= wl_lo) & (wavelengths <= wl_hi)
            band_mask &= ~in_region
        n_excluded = n_bands - np.sum(band_mask)
        if n_excluded > 0:
            print(f"  Excluded {n_excluded} water vapor bands")

    # Apply spectral_opts filtering
    if spectral_opts:
        mode = spectral_opts.get('mode', 'all')

        if mode == 'ranges' and 'ranges' in spectral_opts:
            range_mask = np.zeros(n_bands, dtype=bool)
            for rng in spectral_opts['ranges']:
                lo, hi = rng
                range_mask |= (wavelengths >= lo) & (wavelengths <= hi)
            band_mask &= range_mask
            print(f"  Applied wavelength range filter: "
                  f"{np.sum(band_mask)} bands remaining")

        elif mode == 'nth' and 'nth' in spectral_opts:
            nth = int(spectral_opts['nth'])
            nth_mask = np.zeros(n_bands, dtype=bool)
            valid_indices = np.where(band_mask)[0]
            nth_mask[valid_indices[::nth]] = True
            band_mask = nth_mask
            print(f"  Every-{nth}-band decimation: "
                  f"{np.sum(band_mask)} bands remaining")

    # Apply mask
    band_indices = np.where(band_mask)[0]
    if len(band_indices) < n_bands:
        rad = rad[band_indices, :, :]
        wavelengths = wavelengths[band_indices]
        fwhm_arr = fwhm_arr[band_indices]
        print(f"  Final band count: {len(band_indices)}")

    # --- Bbox clipping of GLT ---
    if bbox is not None:
        glt_line, glt_sample, easting, northing, transform = \
            _clip_glt_to_bbox(glt_line, glt_sample, easting, northing,
                              transform, crs, bbox)

    return {
        'radiance': rad,
        'wavelengths': wavelengths,
        'fwhm': fwhm_arr,
        'glt_line': glt_line,
        'glt_sample': glt_sample,
        'transform': transform,
        'crs': crs,
        'sensor_name': sensor_name,
    }


def _clip_glt_to_bbox(glt_line, glt_sample, easting, northing,
                       transform, crs, bbox):
    """
    Clip GLT arrays to a bounding box.
    bbox is (min_lon, min_lat, max_lon, max_lat) in geographic coords.
    Converts to UTM if needed.
    """
    try:
        from pyproj import Transformer

        # Convert geographic bbox to UTM
        epsg = crs.to_epsg()
        if epsg and epsg != 4326:
            transformer = Transformer.from_crs(4326, epsg, always_xy=True)
            min_e, min_n = transformer.transform(bbox[0], bbox[1])
            max_e, max_n = transformer.transform(bbox[2], bbox[3])
        else:
            min_e, min_n = bbox[0], bbox[1]
            max_e, max_n = bbox[2], bbox[3]

        # Find easting/northing indices within bbox
        e_mask = (easting >= min_e) & (easting <= max_e)
        n_mask = (northing >= min(min_n, max_n)) & \
                 (northing <= max(min_n, max_n))

        if not np.any(e_mask) or not np.any(n_mask):
            print("  WARNING: Bbox does not overlap with data extent")
            return glt_line, glt_sample, easting, northing, transform

        e_idx = np.where(e_mask)[0]
        n_idx = np.where(n_mask)[0]

        e_start, e_end = e_idx[0], e_idx[-1] + 1
        n_start, n_end = n_idx[0], n_idx[-1] + 1

        glt_line = glt_line[n_start:n_end, e_start:e_end]
        glt_sample = glt_sample[n_start:n_end, e_start:e_end]
        easting = easting[e_start:e_end]
        northing = northing[n_start:n_end]

        # Update transform origin
        dx = abs(transform.a)
        dy = transform.e  # negative
        new_ulx = float(easting[0])
        new_uly = float(northing[0])
        transform = Affine(dx, 0.0, new_ulx, 0.0, dy, new_uly)

        print(f"  Bbox clip: {glt_line.shape[1]}×{glt_line.shape[0]} pixels")

    except ImportError:
        print("  WARNING: pyproj not available, bbox clipping skipped")

    return glt_line, glt_sample, easting, northing, transform


# ============================================================================
# GLT Orthorectification
# ============================================================================

def apply_glt(data, glt_line, glt_sample):
    """
    Apply GLT orthorectification to convert swath data to map grid.

    Parameters:
        data: 3D array (bands, swath_lines, swath_samples)
        glt_line: 2D array of source line indices
        glt_sample: 2D array of source sample indices

    Returns:
        3D array (bands, out_lines, out_samples)
    """
    n_bands = data.shape[0]
    out_lines, out_samples = glt_line.shape

    # GLT uses 1-based indexing; 0 or negative = fill
    valid = (glt_line > 0) & (glt_sample > 0)

    # Convert to 0-based
    src_lines = np.clip(glt_line.astype(np.int32) - 1, 0,
                        data.shape[1] - 1)
    src_samples = np.clip(glt_sample.astype(np.int32) - 1, 0,
                          data.shape[2] - 1)

    output = np.full((n_bands, out_lines, out_samples), np.nan,
                     dtype=np.float32)

    for b in range(n_bands):
        band_data = data[b]
        ortho_band = band_data[src_lines, src_samples]
        ortho_band[~valid] = np.nan
        output[b] = ortho_band

    return output


# ============================================================================
# Spatial Resampling
# ============================================================================

def resample_spatial(data, native_res, target_res, transform):
    """
    Resample by block averaging to a coarser resolution.

    Parameters:
        data: 3D array (bands, lines, samples)
        native_res: native pixel size (m)
        target_res: target pixel size (m)
        transform: source Affine transform

    Returns:
        (resampled_data, new_transform)
    """
    factor = int(round(target_res / native_res))
    if factor <= 1:
        return data, transform

    n_bands, lines, samples = data.shape
    new_lines = lines // factor
    new_samples = samples // factor

    print(f"  Resampling: {native_res:.1f}m -> {target_res:.1f}m "
          f"(factor={factor}, {new_lines}×{new_samples} pixels)")

    # Trim to exact multiple
    trimmed = data[:, :new_lines * factor, :new_samples * factor]

    # Block average with NaN handling
    reshaped = trimmed.reshape(n_bands, new_lines, factor,
                               new_samples, factor)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        result = np.nanmean(reshaped, axis=(2, 4)).astype(np.float32)

    # Update transform
    new_transform = Affine(
        transform.a * factor, transform.b * factor, transform.c,
        transform.d * factor, transform.e * factor, transform.f
    )

    return result, new_transform


# ============================================================================
# Output writers
# ============================================================================

def write_output(data_3d, output_path, crs, transform,
                 band_names=None, wavelengths=None, fwhm=None,
                 description='', write_geotiff=True, write_envi=True):
    """
    Write multiband output as GeoTIFF and/or ENVI.
    Input: (bands, height, width) — already in bands-first order.
    """
    output_files = []
    num_bands, height, width = data_3d.shape

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
              f"({num_bands} bands, {width}×{height})")
        with rasterio.open(gtiff_path, 'w', **profile) as dst:
            for i in range(num_bands):
                dst.write(data_3d[i], i + 1)
                if band_names and i < len(band_names):
                    dst.set_band_description(i + 1, band_names[i])
        output_files.append(gtiff_path)

    if write_envi:
        envi_dat = output_path + '.dat'
        envi_hdr = output_path + '.hdr'
        print(f"    Writing ENVI: {os.path.basename(envi_dat)} "
              f"({num_bands} bands, {width}×{height})")

        # Replace NaN with -9999 for ENVI
        envi_data = data_3d.copy()
        envi_data[np.isnan(envi_data)] = -9999

        envi_data.astype('<f4').tofile(envi_dat)
        _build_envi_header(
            envi_hdr, height, width, num_bands, crs, transform,
            wavelengths=wavelengths, fwhm=fwhm, band_names=band_names,
            description=description
        )
        output_files.append(envi_dat)
        output_files.append(envi_hdr)

    return output_files


def _build_envi_header(hdr_path, height, width, num_bands, crs, transform,
                       wavelengths=None, fwhm=None, band_names=None,
                       description=''):
    """Write ENVI header with wavelength metadata."""
    # Build map info
    epsg = crs.to_epsg() if crs else None
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
        map_info = (f'UTM, 1, 1, {ulx}, {uly}, {dx}, {dy}, '
                    f'{zone}, {hemi}, WGS-84, units=Meters')
    else:
        map_info = f'Geographic Lat/Lon, 1, 1, {ulx}, {uly}, {dx}, {dy}'

    crs_wkt = crs.to_wkt() if crs else ''

    lines = [
        'ENVI',
        f'description = {{{description}}}',
        f'samples = {width}',
        f'lines   = {height}',
        f'bands   = {num_bands}',
        'header offset = 0',
        'file type = ENVI Standard',
        'data type = 4',
        'interleave = bsq',
        'byte order = 0',
        f'map info = {{{map_info}}}',
    ]

    if crs_wkt:
        lines.append(f'coordinate system string = {{{crs_wkt}}}')

    if wavelengths is not None and len(wavelengths) > 0:
        wl_str = ', '.join(f'{w:.4f}' for w in wavelengths)
        lines.append(f'wavelength = {{{wl_str}}}')
        lines.append('wavelength units = Nanometers')

    if fwhm is not None and len(fwhm) > 0:
        fwhm_str = ', '.join(f'{f:.4f}' for f in fwhm)
        lines.append(f'fwhm = {{{fwhm_str}}}')

    if band_names:
        bn_str = ', '.join(band_names)
        lines.append(f'band names = {{{bn_str}}}')

    lines.append('data ignore value = -9999')

    with open(hdr_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


# ============================================================================
# Main processing
# ============================================================================

def process_granule(filepath, granule_id, output_dir, exclude_water=True,
                    spectral_opts=None, bbox=None, resolution_m=None,
                    write_geotiff=True, write_envi=True, overwrite=True):
    """Process a single AVIRIS-3/5 L1B granule."""
    all_outputs = []

    data = read_av_data(filepath, exclude_water=exclude_water,
                        spectral_opts=spectral_opts, bbox=bbox)

    sensor_name = data['sensor_name']
    wl = data['wavelengths']
    radiance = data['radiance']

    print(f"\n  Applying GLT orthorectification...")
    ortho = apply_glt(radiance, data['glt_line'], data['glt_sample'])
    del radiance

    out_transform = data['transform']

    # Apply spatial resampling if requested
    native_res = abs(data['transform'].a)
    if resolution_m is not None and resolution_m > native_res:
        ortho, out_transform = resample_spatial(
            ortho, native_res, resolution_m, data['transform'])

    os.makedirs(output_dir, exist_ok=True)

    band_names = [f'{w:.1f} nm' for w in wl]

    # Use sensor prefix in output filename
    prefix = sensor_name.replace('-', '')  # AV3 -> AVIRIS3, AV5 -> AVIRIS5
    if sensor_name == 'AVIRIS-3':
        prefix = 'AV3'
    elif sensor_name == 'AVIRIS-5':
        prefix = 'AV5'

    rad_base = os.path.join(output_dir,
                            f"{prefix}_L1B_Radiance_{granule_id}")

    # Check overwrite
    tif_path = rad_base + '.tif'
    dat_path = rad_base + '.dat'
    if not overwrite:
        tif_exists = os.path.isfile(tif_path) if write_geotiff else True
        dat_exists = os.path.isfile(dat_path) if write_envi else True
        if tif_exists and dat_exists:
            print(f"  Output exists, skipping (overwrite=False)")
            return all_outputs

    outputs = write_output(
        ortho, rad_base, data['crs'], out_transform,
        band_names=band_names,
        wavelengths=wl,
        fwhm=data['fwhm'],
        description=f'{sensor_name} L1B Calibrated Radiance ({granule_id})',
        write_geotiff=write_geotiff, write_envi=write_envi
    )
    all_outputs.extend(outputs)

    return all_outputs


def process_directory(input_dir, output_dir=None, exclude_water=True,
                      spectral_opts=None, bbox=None, resolution_m=None,
                      write_geotiff=True, write_envi=True, overwrite=True):
    """
    Process all AVIRIS-3/5 L1B files in a directory.

    Parameters:
        input_dir: directory containing AVIRIS-3/5 NetCDF files
        output_dir: output directory (default: input_dir/AVIRIS35_output)
        exclude_water: auto-exclude water vapor absorption bands
        spectral_opts: dict with spectral filtering options
        bbox: bounding box (min_lon, min_lat, max_lon, max_lat)
        resolution_m: target pixel size for block-average resampling (m)
        write_geotiff: output GeoTIFF
        write_envi: output ENVI
        overwrite: overwrite existing output files
    """
    if output_dir is None:
        output_dir = os.path.join(input_dir, 'AVIRIS35_output')

    if not write_geotiff and not write_envi:
        write_geotiff = True

    print(f"AVIRIS-3/5 L1B Radiance Converter")
    print(f"{'='*70}")
    print(f"Input:  {input_dir}")
    print(f"Output: {output_dir}")
    fmt = []
    if write_geotiff:
        fmt.append('GeoTIFF')
    if write_envi:
        fmt.append('ENVI')
    print(f"Format: {'/'.join(fmt)}")

    if exclude_water:
        print(f"Water vapor exclusion: ON")
    if resolution_m:
        print(f"Target resolution: {resolution_m} m")

    files = discover_files(input_dir)

    if not files:
        print("\nNo AVIRIS-3 or AVIRIS-5 L1B files found.")
        print("Expected patterns:")
        print("  AV3{YYYYMMDD}t{HHMMSS}*L1B*RDN*.nc  (AVIRIS-3)")
        print("  AV5{YYYYMMDD}t{HHMMSS}*L1B*RDN*.nc  (AVIRIS-5)")
        return

    # Count by sensor
    av3_count = sum(1 for _, gid in files if gid.upper().startswith('AV3'))
    av5_count = sum(1 for _, gid in files if gid.upper().startswith('AV5'))
    print(f"\nFound {len(files)} file(s)", end='')
    parts = []
    if av3_count:
        parts.append(f"{av3_count} AVIRIS-3")
    if av5_count:
        parts.append(f"{av5_count} AVIRIS-5")
    if parts:
        print(f": {', '.join(parts)}")
    else:
        print()

    for i, (filepath, granule_id) in enumerate(files, 1):
        print(f"\n[{i}/{len(files)}] {granule_id}")
        try:
            process_granule(
                filepath, granule_id, output_dir,
                exclude_water=exclude_water,
                spectral_opts=spectral_opts,
                bbox=bbox,
                resolution_m=resolution_m,
                write_geotiff=write_geotiff,
                write_envi=write_envi,
                overwrite=overwrite,
            )
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()

    print(f"\n{'='*70}")
    print(f"Processing complete. Output: {output_dir}")


# ============================================================================
# CLI
# ============================================================================

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='AVIRIS-3/5 L1B Converter — '
                    'Convert NetCDF radiance data to GeoTIFF/ENVI')
    parser.add_argument('input_dir',
                        help='Input directory with AVIRIS-3/5 .nc files')
    parser.add_argument('-o', '--output', help='Output directory')
    parser.add_argument('--geotiff', action='store_true', default=True,
                        help='Write GeoTIFF (default)')
    parser.add_argument('--envi', action='store_true',
                        help='Write ENVI output')
    parser.add_argument('--both', action='store_true',
                        help='Write both GeoTIFF and ENVI')
    parser.add_argument('--no-water-exclude', action='store_true',
                        help='Keep water vapor absorption bands')
    parser.add_argument('--resolution', type=float,
                        help='Target resolution for resampling (m)')
    parser.add_argument('--no-overwrite', action='store_true',
                        help='Skip existing output files')

    args = parser.parse_args()

    do_geotiff = True
    do_envi = args.envi
    if args.both:
        do_geotiff = True
        do_envi = True

    process_directory(
        args.input_dir,
        output_dir=args.output,
        exclude_water=not args.no_water_exclude,
        resolution_m=args.resolution,
        write_geotiff=do_geotiff,
        write_envi=do_envi,
        overwrite=not args.no_overwrite,
    )
