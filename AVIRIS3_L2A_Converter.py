"""
AVIRIS-3/5 L2A Surface Reflectance Converter
===============================================
Reads AVIRIS-3/5 L2A orthorectified NetCDF files (from EarthData/ORNL DAAC)
and outputs multiband GeoTIFF and/or ENVI files.

The L2A files are already orthorectified on a UTM grid — no GLT or
reprojection is needed.  Data is read from NetCDF groups:
  - reflectance group: reflectance (nbands, nrows, ncols), wavelength, fwhm
  - aerosol_optical_thickness group: single-band AOT
  - water_vapor group: single-band water vapor (cm)

Output files:
  - Surface Reflectance (multiband, water vapor bands optionally excluded)
  - Ancillary (AOT + Water Vapor, 2 bands) [optional]

Filename pattern:
  AV3{YYYYMMDD}t{HHMMSS}_{scene}_L2A_OE_{hash}_RFL_ORT.nc  (AVIRIS-3)
  AV5{YYYYMMDD}t{HHMMSS}_{scene}_L2A_OE_{hash}_RFL_ORT.nc  (AVIRIS-5)

Dependencies: h5py (or netCDF4), numpy, rasterio
"""

import os
import re
import sys
import numpy as np

try:
    import rasterio
    from rasterio.crs import CRS
    from rasterio.transform import Affine
except ImportError:
    print("ERROR: rasterio is required.")
    sys.exit(1)

# Try netCDF4 first, fall back to h5py
_NC_BACKEND = None
try:
    import netCDF4
    _NC_BACKEND = 'netcdf4'
except ImportError:
    try:
        import h5py
        _NC_BACKEND = 'h5py'
    except ImportError:
        print("ERROR: netCDF4 or h5py is required.")
        sys.exit(1)


# ============================================================================
# Constants
# ============================================================================

# Filename patterns for L2A reflectance files
AV_L2A_RFL_PATTERN = re.compile(
    r'^(AV[35]\d{8}t\d{6}).*L2A.*_RFL_ORT\.nc$',
    re.IGNORECASE
)

# Also match uncertainty files
AV_L2A_UNC_PATTERN = re.compile(
    r'^(AV[35]\d{8}t\d{6}).*L2A.*_UNC_ORT\.nc$',
    re.IGNORECASE
)

# Water vapor absorption regions to auto-exclude (nm)
WATER_VAPOR_REGIONS = [
    (1334, 1440),   # 1.4 um H2O (extended to catch edge channels)
    (1788, 1980),   # 1.9 um H2O
]


# ============================================================================
# NetCDF reading abstraction
# ============================================================================

def nc_open(filepath):
    """Open a NetCDF/HDF5 file."""
    if _NC_BACKEND == 'netcdf4':
        return netCDF4.Dataset(filepath, 'r')
    else:
        return h5py.File(filepath, 'r')


def nc_close(ds):
    """Close the dataset."""
    ds.close()


def nc_get_group(ds, group_name):
    """Get a group from the dataset."""
    if _NC_BACKEND == 'netcdf4':
        return ds.groups[group_name]
    else:
        return ds[group_name]


def nc_read_var(group, varname):
    """Read a variable as numpy array (raw, no auto-scaling)."""
    if _NC_BACKEND == 'netcdf4':
        var = group.variables[varname]
        var.set_auto_maskandscale(False)
        return np.array(var[:])
    else:
        return np.array(group[varname])


def nc_get_attr(ds_or_group, varname_or_attr, attr=None, default=None):
    """Get an attribute. If attr is None, get a global/group attribute."""
    try:
        if attr is not None:
            # Variable attribute
            if _NC_BACKEND == 'netcdf4':
                return ds_or_group.variables[varname_or_attr].getncattr(attr)
            else:
                return ds_or_group[varname_or_attr].attrs[attr]
        else:
            # Global/group attribute
            if _NC_BACKEND == 'netcdf4':
                return ds_or_group.getncattr(varname_or_attr)
            else:
                return ds_or_group.attrs[varname_or_attr]
    except (KeyError, AttributeError):
        return default


def nc_get_global_attr(ds, attrname, default=None):
    """Get a global attribute from root dataset."""
    try:
        if _NC_BACKEND == 'netcdf4':
            return ds.getncattr(attrname)
        else:
            return ds.attrs[attrname]
    except (KeyError, AttributeError):
        return default


# ============================================================================
# File discovery
# ============================================================================

def discover_files(input_dir):
    """
    Scan for AVIRIS-3/5 L2A RFL_ORT files.
    Returns list of (rfl_path, unc_path_or_None, scene_id).
    """
    rfl_files = {}
    unc_files = {}

    def scan_dir(d):
        try:
            entries = os.listdir(d)
        except OSError:
            return
        for f in entries:
            m = AV_L2A_RFL_PATTERN.match(f)
            if m:
                gid = m.group(1)
                # Include scene number in ID if present
                scene_match = re.search(r'_(\d{3})_L2A', f)
                scene = scene_match.group(1) if scene_match else '001'
                key = f"{gid}_{scene}"
                rfl_files[key] = os.path.join(d, f)
            m = AV_L2A_UNC_PATTERN.match(f)
            if m:
                gid = m.group(1)
                scene_match = re.search(r'_(\d{3})_L2A', f)
                scene = scene_match.group(1) if scene_match else '001'
                key = f"{gid}_{scene}"
                unc_files[key] = os.path.join(d, f)

    scan_dir(input_dir)
    for entry in os.listdir(input_dir):
        subdir = os.path.join(input_dir, entry)
        if os.path.isdir(subdir):
            scan_dir(subdir)

    # Pair RFL with UNC where available
    results = []
    for key in sorted(rfl_files.keys()):
        rfl_path = rfl_files[key]
        unc_path = unc_files.get(key)
        results.append((rfl_path, unc_path, key))

    return results


# ============================================================================
# CRS and transform extraction
# ============================================================================

def extract_crs_transform(ds):
    """
    Extract CRS and affine transform from the transverse_mercator variable
    or global attributes.
    """
    crs = None
    transform = None

    # Try reading from transverse_mercator variable
    try:
        if _NC_BACKEND == 'netcdf4':
            tm_var = ds.variables['transverse_mercator']
            wkt = tm_var.getncattr('crs_wkt') if 'crs_wkt' in tm_var.ncattrs() else None
            if wkt is None:
                wkt = tm_var.getncattr('spatial_ref') if 'spatial_ref' in tm_var.ncattrs() else None
            geo_str = tm_var.getncattr('GeoTransform') if 'GeoTransform' in tm_var.ncattrs() else None
        else:
            tm_var = ds['transverse_mercator']
            wkt = tm_var.attrs.get('crs_wkt', tm_var.attrs.get('spatial_ref', None))
            geo_str = tm_var.attrs.get('GeoTransform', None)

        if wkt is not None:
            if isinstance(wkt, bytes):
                wkt = wkt.decode('utf-8')
            crs = CRS.from_wkt(wkt)

        if geo_str is not None:
            if isinstance(geo_str, bytes):
                geo_str = geo_str.decode('utf-8')
            parts = [float(x) for x in str(geo_str).split()]
            # GeoTransform = ul_x, pixel_width, rot_x, ul_y, rot_y, pixel_height
            transform = Affine(parts[1], parts[2], parts[0],
                               parts[4], parts[5], parts[3])
    except (KeyError, IndexError):
        pass

    # Fallback: build from easting/northing coordinate arrays
    if transform is None:
        try:
            if _NC_BACKEND == 'netcdf4':
                easting = np.array(ds.variables['easting'][:])
                northing = np.array(ds.variables['northing'][:])
            else:
                easting = np.array(ds['easting'])
                northing = np.array(ds['northing'])

            pixel_x = float(easting[1] - easting[0])
            pixel_y = float(northing[1] - northing[0])
            # Coordinates are cell centers, transform needs upper-left corner
            ul_x = float(easting[0]) - pixel_x / 2
            ul_y = float(northing[0]) - pixel_y / 2
            transform = Affine(pixel_x, 0, ul_x, 0, pixel_y, ul_y)
        except (KeyError, IndexError):
            pass

    return crs, transform


# ============================================================================
# ENVI header writing
# ============================================================================

def build_envi_header(hdr_path, height, width, num_bands, dtype, crs,
                      transform, wavelengths=None, fwhm=None,
                      band_names=None, description='', interleave='bsq'):
    """Write a custom ENVI header file with wavelength metadata."""
    dtype_map = {
        'uint8': 1, 'int16': 2, 'int32': 3, 'float32': 4,
        'float64': 5, 'uint16': 12, 'uint32': 13,
    }
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

        if transform is not None and crs is not None:
            x_size = abs(transform.a)
            y_size = abs(transform.e)
            ul_x = transform.c
            ul_y = transform.f

            if crs.is_geographic:
                f.write(f'map info = {{Geographic Lat/Lon, 1, 1, {ul_x}, '
                        f'{ul_y}, {x_size}, {y_size}, WGS-84}}\n')
            elif crs.is_projected:
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
            wl_str = ', '.join(f'{w:.2f}' for w in wavelengths)
            f.write(f'wavelength = {{\n  {wl_str}}}\n')

        if fwhm is not None:
            fwhm_str = ', '.join(f'{w:.2f}' for w in fwhm)
            f.write(f'fwhm = {{\n  {fwhm_str}}}\n')

        if band_names is not None:
            bn_str = ',\n  '.join(band_names)
            f.write(f'band names = {{\n  {bn_str}}}\n')

        f.write('data ignore value = -9999.0\n')


# ============================================================================
# Output writing
# ============================================================================

def write_multiband(stack, output_dir, basename, crs, transform,
                    wavelengths=None, fwhm=None, band_names=None,
                    description='', write_geotiff=True, write_envi=True):
    """Write 3D array (bands, rows, cols) as GeoTIFF and/or ENVI."""
    output_files = []
    num_bands, height, width = stack.shape

    if write_geotiff:
        gtiff_dir = os.path.join(output_dir, 'GeoTIFF')
        os.makedirs(gtiff_dir, exist_ok=True)
        gtiff_path = os.path.join(gtiff_dir, basename + '.tif')
        profile = {
            'driver': 'GTiff', 'dtype': 'float32',
            'width': width, 'height': height, 'count': num_bands,
            'crs': crs, 'transform': transform, 'nodata': -9999.0,
            'compress': 'deflate', 'predictor': 2, 'zlevel': 6,
            'tiled': True, 'blockxsize': 256, 'blockysize': 256,
        }
        print(f"    Writing GeoTIFF: {basename}.tif")
        with rasterio.open(gtiff_path, 'w', **profile) as dst:
            for i in range(num_bands):
                dst.write(stack[i], i + 1)
                if band_names and i < len(band_names):
                    dst.set_band_description(i + 1, band_names[i])
        output_files.append(gtiff_path)

    if write_envi:
        envi_dir = os.path.join(output_dir, 'ENVI')
        os.makedirs(envi_dir, exist_ok=True)
        envi_dat = os.path.join(envi_dir, basename + '.dat')
        envi_hdr = os.path.join(envi_dir, basename + '.hdr')
        print(f"    Writing ENVI: {basename}.dat")
        stack.tofile(envi_dat)
        build_envi_header(
            envi_hdr, height=height, width=width, num_bands=num_bands,
            dtype=np.float32, crs=crs, transform=transform,
            wavelengths=wavelengths, fwhm=fwhm, band_names=band_names,
            description=description, interleave='bsq')
        output_files.append(envi_dat)
        output_files.append(envi_hdr)

    return output_files


# ============================================================================
# Band filtering
# ============================================================================

def filter_water_vapor(wavelengths, exclude_regions=None):
    """Return boolean mask of bands to keep (True = keep)."""
    if exclude_regions is None:
        exclude_regions = WATER_VAPOR_REGIONS
    keep = np.ones(len(wavelengths), dtype=bool)
    for wmin, wmax in exclude_regions:
        keep &= ~((wavelengths >= wmin) & (wavelengths <= wmax))
    return keep


def apply_spectral_filter(wavelengths, spectral_opts):
    """Apply spectral subsetting options. Returns boolean mask."""
    keep = np.ones(len(wavelengths), dtype=bool)
    if spectral_opts is None:
        return keep

    wl_min = spectral_opts.get('wl_min')
    wl_max = spectral_opts.get('wl_max')
    if wl_min is not None:
        keep &= (wavelengths >= wl_min)
    if wl_max is not None:
        keep &= (wavelengths <= wl_max)

    step = spectral_opts.get('band_step')
    if step and step > 1:
        step_mask = np.zeros(len(wavelengths), dtype=bool)
        step_mask[::step] = True
        keep &= step_mask

    return keep


# ============================================================================
# Main processing
# ============================================================================

def process_file(rfl_path, unc_path, scene_id, output_dir,
                 exclude_water=True, spectral_opts=None, bbox=None,
                 write_geotiff=True, write_envi=True, overwrite=True):
    """Process a single AVIRIS-3/5 L2A file.
    bbox: (min_lon, min_lat, max_lon, max_lat) for spatial subsetting.
    """
    all_outputs = []
    fname = os.path.basename(rfl_path)

    # Detect sensor
    sensor = 'AV3' if scene_id.upper().startswith('AV3') else 'AV5'
    sensor_label = 'AVIRIS-3' if sensor == 'AV3' else 'AVIRIS-5'

    print(f"\n  Sensor: {sensor_label}")
    print(f"  File: {fname}")

    # Open and extract CRS/transform
    ds = nc_open(rfl_path)
    crs, transform = extract_crs_transform(ds)

    if crs is None or transform is None:
        print("  ERROR: Could not extract CRS or transform.")
        nc_close(ds)
        return all_outputs

    pixel_size = abs(transform.a)
    print(f"  CRS: {crs.to_epsg() or 'custom UTM'}")
    print(f"  Pixel size: {pixel_size:.1f} m")

    # Read reflectance group
    try:
        rfl_group = nc_get_group(ds, 'reflectance')
    except (KeyError, ValueError):
        print("  ERROR: 'reflectance' group not found in file.")
        nc_close(ds)
        return all_outputs

    wavelengths = nc_read_var(rfl_group, 'wavelength')
    fwhm_arr = nc_read_var(rfl_group, 'fwhm')
    rfl_data = nc_read_var(rfl_group, 'reflectance')  # (nbands, nrows, ncols)

    fill_value = -9999.0
    try:
        if _NC_BACKEND == 'netcdf4':
            fill_value = float(rfl_group.variables['reflectance'].getncattr('_FillValue'))
        else:
            fill_value = float(rfl_group['reflectance'].attrs['_FillValue'])
    except (KeyError, AttributeError):
        pass

    nbands, nrows, ncols = rfl_data.shape
    print(f"  Dimensions: {nbands} bands x {nrows} rows x {ncols} cols")
    print(f"  Wavelength range: {wavelengths[0]:.1f} - {wavelengths[-1]:.1f} nm")

    # Read easting/northing for bbox cropping
    try:
        if _NC_BACKEND == 'netcdf4':
            easting = np.array(ds.variables['easting'][:])
            northing = np.array(ds.variables['northing'][:])
        else:
            easting = np.array(ds['easting'])
            northing = np.array(ds['northing'])
    except KeyError:
        easting = None
        northing = None

    # Apply spatial subset (bbox in lon/lat -> UTM -> row/col crop)
    row_slice = slice(None)
    col_slice = slice(None)

    if bbox is not None and easting is not None and northing is not None:
        min_lon, min_lat, max_lon, max_lat = bbox
        try:
            from pyproj import Transformer
            transformer = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
            min_x, min_y = transformer.transform(min_lon, min_lat)
            max_x, max_y = transformer.transform(max_lon, max_lat)

            # Find column range (easting)
            col_mask = (easting >= min(min_x, max_x)) & \
                       (easting <= max(min_x, max_x))
            if col_mask.any():
                col_indices = np.where(col_mask)[0]
                col_slice = slice(int(col_indices[0]), int(col_indices[-1]) + 1)

            # Find row range (northing)
            row_mask = (northing >= min(min_y, max_y)) & \
                       (northing <= max(min_y, max_y))
            if row_mask.any():
                row_indices = np.where(row_mask)[0]
                row_slice = slice(int(row_indices[0]), int(row_indices[-1]) + 1)

            print(f"  Spatial subset: rows [{row_slice.start}:{row_slice.stop}], "
                  f"cols [{col_slice.start}:{col_slice.stop}]")
        except ImportError:
            print("  WARNING: pyproj not available, skipping bbox subset.")
        except Exception as e:
            print(f"  WARNING: Bbox conversion failed ({e}), skipping subset.")

    # Crop data and coordinates
    rfl_data = rfl_data[:, row_slice, col_slice]
    nbands, nrows, ncols = rfl_data.shape

    # Update transform for cropped region
    if row_slice != slice(None) or col_slice != slice(None):
        crop_easting = easting[col_slice]
        crop_northing = northing[row_slice]
        pixel_x = float(easting[1] - easting[0]) if len(easting) > 1 else pixel_size
        pixel_y = float(northing[1] - northing[0]) if len(northing) > 1 else -pixel_size
        ul_x = float(crop_easting[0]) - abs(pixel_x) / 2
        ul_y = float(crop_northing[0]) - pixel_y / 2
        transform = Affine(pixel_x, 0, ul_x, 0, pixel_y, ul_y)
        print(f"  Cropped dimensions: {nbands} bands x {nrows} rows x {ncols} cols")

    # Apply band filtering
    keep = np.ones(nbands, dtype=bool)

    if exclude_water:
        water_mask = filter_water_vapor(wavelengths)
        n_excluded = np.sum(~water_mask)
        if n_excluded > 0:
            print(f"  Excluding {n_excluded} water vapor bands")
        keep &= water_mask

    if spectral_opts:
        spec_mask = apply_spectral_filter(wavelengths, spectral_opts)
        keep &= spec_mask

    # Apply mask
    sel_indices = np.where(keep)[0]
    n_out = len(sel_indices)
    print(f"  Output bands: {n_out}")

    if n_out == 0:
        print("  WARNING: No bands remaining after filtering.")
        nc_close(ds)
        return all_outputs

    sel_wavelengths = wavelengths[sel_indices]
    sel_fwhm = fwhm_arr[sel_indices]
    sel_data = rfl_data[sel_indices, :, :]

    # Mask fill values
    sel_data = np.float32(sel_data)
    sel_data[sel_data == fill_value] = np.nan

    # Build band names
    band_names = [f'{sel_wavelengths[i]:.1f} nm' for i in range(n_out)]

    # Build output basename
    out_basename = f"{sensor}_{scene_id}_L2A_SurfRefl"

    # Check overwrite
    gtiff_check = os.path.join(output_dir, 'GeoTIFF', out_basename + '.tif')
    if not overwrite and os.path.isfile(gtiff_check):
        print("  Output exists, skipping (overwrite=False)")
        nc_close(ds)
        return all_outputs

    # Write reflectance
    outputs = write_multiband(
        sel_data, output_dir, out_basename, crs, transform,
        wavelengths=sel_wavelengths, fwhm=sel_fwhm,
        band_names=band_names,
        description=f'{sensor_label} L2A Surface Reflectance ({scene_id})',
        write_geotiff=write_geotiff, write_envi=write_envi)
    all_outputs.extend(outputs)

    # Read ancillary data (AOT + Water Vapor)
    anc_data = []
    anc_names = []

    try:
        aot_group = nc_get_group(ds, 'aerosol_optical_thickness')
        aot = nc_read_var(aot_group, 'aerosol_optical_thickness')
        aot = np.float32(aot[row_slice, col_slice])
        aot[aot == fill_value] = np.nan
        anc_data.append(aot)
        anc_names.append('Aerosol Optical Thickness')
    except (KeyError, ValueError):
        pass

    try:
        wv_group = nc_get_group(ds, 'water_vapor')
        wv = nc_read_var(wv_group, 'water_vapor')
        wv = np.float32(wv[row_slice, col_slice])
        wv[wv == fill_value] = np.nan
        anc_data.append(wv)
        anc_names.append('Water Vapor (cm)')
    except (KeyError, ValueError):
        pass

    if anc_data:
        anc_stack = np.stack(anc_data, axis=0)
        anc_basename = f"{sensor}_{scene_id}_L2A_Ancillary"
        outputs = write_multiband(
            anc_stack, output_dir, anc_basename, crs, transform,
            band_names=anc_names,
            description=f'{sensor_label} L2A Ancillary ({scene_id})',
            write_geotiff=write_geotiff, write_envi=write_envi)
        all_outputs.extend(outputs)

    nc_close(ds)
    return all_outputs


def process_directory(input_dir, output_dir=None, exclude_water=True,
                      spectral_opts=None, bbox=None, write_geotiff=True,
                      write_envi=True, overwrite=True):
    """Process all AVIRIS-3/5 L2A files in a directory."""
    if output_dir is None:
        output_dir = os.path.join(input_dir, 'AVIRIS35_L2A_output')

    print(f"AVIRIS-3/5 L2A Surface Reflectance Converter")
    print(f"{'='*70}")
    print(f"Input:  {input_dir}")
    print(f"Output: {output_dir}")

    files = discover_files(input_dir)

    if not files:
        print("\nNo AVIRIS-3/5 L2A files found.")
        print("Expected pattern: AV3/AV5*_L2A_*_RFL_ORT.nc")
        return

    av3_count = sum(1 for _, _, gid in files if gid.upper().startswith('AV3'))
    av5_count = len(files) - av3_count
    parts = []
    if av3_count:
        parts.append(f"{av3_count} AVIRIS-3")
    if av5_count:
        parts.append(f"{av5_count} AVIRIS-5")
    print(f"\nFound {len(files)} file(s): {', '.join(parts)}\n")

    os.makedirs(output_dir, exist_ok=True)
    all_outputs = []

    for i, (rfl_path, unc_path, scene_id) in enumerate(files, 1):
        print(f"{'='*70}")
        print(f"File {i}/{len(files)}: {scene_id}")
        print(f"{'='*70}")

        try:
            outputs = process_file(
                rfl_path, unc_path, scene_id, output_dir,
                exclude_water=exclude_water, spectral_opts=spectral_opts,
                bbox=bbox,
                write_geotiff=write_geotiff, write_envi=write_envi,
                overwrite=overwrite)
            all_outputs.extend(outputs)
        except Exception as e:
            print(f"  ERROR: {e}")
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
        print("AVIRIS-3/5 L2A Surface Reflectance Converter")
        print("=" * 50)
        print(f"\nUsage: python {os.path.basename(__file__)} <input_dir> "
              f"[output_dir] [options]")
        print(f"\nSupported: AVIRIS-3 and AVIRIS-5 L2A ORT NetCDF files")
        print(f"\nOptions:")
        print(f"  --geotiff    Output GeoTIFF only")
        print(f"  --envi       Output ENVI only")
        print(f"  --keep-water Keep water vapor absorption bands")
        print(f"  (default: both formats, water vapor excluded)")
        sys.exit(0)

    args = sys.argv[1:]
    write_geotiff = True
    write_envi = True
    exclude_water = True
    positional = []

    for arg in args:
        if arg.lower() == '--geotiff':
            write_geotiff, write_envi = True, False
        elif arg.lower() == '--envi':
            write_geotiff, write_envi = False, True
        elif arg.lower() == '--keep-water':
            exclude_water = False
        else:
            positional.append(arg)

    input_dir = positional[0]
    output_dir = positional[1] if len(positional) > 1 else None

    if not os.path.isdir(input_dir):
        print(f"Error: {input_dir} does not exist")
        sys.exit(1)

    process_directory(input_dir, output_dir,
                      exclude_water=exclude_water,
                      write_geotiff=write_geotiff, write_envi=write_envi)
