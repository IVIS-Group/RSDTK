"""
ASTER_Converter.py  —  part of the Remote Sensing Data Toolkit (RSDTK)

Converts ASTER products to analysis-ready multiband GeoTIFF and/or ENVI files
with scaling factors applied in a single pass.

Supported V004 GeoTIFF products:
    AST_05   — TIR Surface Emissivity         (÷ 1000  → 0–1)
    AST_07   — Surface Reflectance VNIR/SWIR  (÷ 1000  → 0–1)  [auto-detects 07XT]
    AST_08   — Surface Kinetic Temperature    (÷ 10    → K)
    AST_09T  — TIR Calibrated Radiance        (× per-band UCC  → W/m²/sr/µm)
               Outputs two files: SurfaceRadiance and SkyIrradiance
    AST_L1T  — Precision Terrain Corrected Radiance
               VNIR (bands 1-3N): raw DN  [no scaling]
               SWIR (bands 4-9):  raw DN  [no scaling, optional — may be absent]
               TIR  (bands 10-14): Radiance = (DN - 1) × UCC  [W/m²/sr/µm]
               Outputs three files: VNIR, SWIR (if present), TIR

Supported V003 HDF products (requires pyhdf):
    AST_05   — TIR Surface Emissivity         (× 0.001 → 0–1)
    AST_07   — Surface Reflectance VNIR       (× 0.001 → 0–1)
    AST_08   — Surface Kinetic Temperature    (× 0.1   → K)
    AST_09T  — TIR Calibrated Radiance        (× per-band UCC  → W/m²/sr/µm)

    V003 geolocation: .aux.xml (AST_05, AST_09T) or .met corner coordinates.

Input layout (both supported automatically):
    • Subfolder-per-granule  (standard NASA LP DAAC download)
    • Flat folder            (all granule files in one directory)

Output:
    <output_dir>/GeoTIFF/<granule_id>_<product>.tif
    <output_dir>/ENVI/<granule_id>_<product>.dat + .hdr

Usage (CLI):
    python ASTER_Converter.py <input_dir> [output_dir] --product AST_05 [--geotiff|--envi]
    python ASTER_Converter.py <input_dir> --product AST_09T --v003

Usage (RSDTK):
    from ASTER_Converter import process_directory, process_directory_v003
    process_directory(input_dir, output_dir, aster_product='AST_05',
                      write_geotiff=True, write_envi=True)
    process_directory_v003(input_dir, output_dir, aster_product='AST_05',
                           write_geotiff=True, write_envi=True)
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
    print("ERROR: rasterio is required.  conda install -c conda-forge rasterio")
    sys.exit(1)

# V003 HDF4 support (optional — only needed for .hdf files)
try:
    from pyhdf.SD import SD, SDC
    _HAS_PYHDF = True
except ImportError:
    _HAS_PYHDF = False


# =============================================================================
# PRODUCT DEFINITIONS
# =============================================================================

# Central wavelengths (µm) for each ASTER band number
ASTER_WAVELENGTHS = {
    1:  0.556,
    2:  0.661,
    3:  0.807,
    4:  1.656,
    5:  2.167,
    6:  2.209,
    7:  2.262,
    8:  2.336,
    9:  2.400,
    10: 8.291,
    11: 8.634,
    12: 9.075,
    13: 10.657,
    14: 11.318,
}

# DN-to-radiance unit conversion coefficients for AST_09T and AST_L1T TIR
# Radiance = (DN - 1) × UCC   [W/(m² · sr · µm)]
TIR_UCC = {
    10: 0.006882,
    11: 0.006780,
    12: 0.006590,
    13: 0.005693,
    14: 0.005225,
}

ASTER_PRODUCTS = {
    "AST_05": {
        "name": "TIR Surface Emissivity",
        "science_bands": {
            "Band_10": {"aster_band": 10, "description": "Surface Emissivity Band 10"},
            "Band_11": {"aster_band": 11, "description": "Surface Emissivity Band 11"},
            "Band_12": {"aster_band": 12, "description": "Surface Emissivity Band 12"},
            "Band_13": {"aster_band": 13, "description": "Surface Emissivity Band 13"},
            "Band_14": {"aster_band": 14, "description": "Surface Emissivity Band 14"},
        },
        "scaling": {"method": "divide", "factor": 1000},
        "valid_range": (0.0, 1.0),
        "output_suffix": "",
    },
    "AST_07X": {
        "name": "Surface Reflectance (VNIR/SWIR)",
        "science_bands": {
            "Band_01": {"aster_band": 1,    "description": "VNIR Reflectance Band 1"},
            "Band_02": {"aster_band": 2,    "description": "VNIR Reflectance Band 2"},
            "Band_3N": {"aster_band": 3,    "description": "VNIR Reflectance Band 3N"},
            "Band_04": {"aster_band": 4,    "description": "SWIR Reflectance Band 4",  "optional": True},
            "Band_05": {"aster_band": 5,    "description": "SWIR Reflectance Band 5",  "optional": True},
            "Band_06": {"aster_band": 6,    "description": "SWIR Reflectance Band 6",  "optional": True},
            "Band_07": {"aster_band": 7,    "description": "SWIR Reflectance Band 7",  "optional": True},
            "Band_08": {"aster_band": 8,    "description": "SWIR Reflectance Band 8",  "optional": True},
            "Band_09": {"aster_band": 9,    "description": "SWIR Reflectance Band 9",  "optional": True},
        },
        "scaling": {"method": "divide", "factor": 1000},
        "valid_range": (0.0, 1.0),
        "output_suffix": "",
    },
    "AST_08": {
        "name": "Surface Kinetic Temperature",
        "science_bands": {
            "Temperature": {"aster_band": None, "description": "Surface Kinetic Temperature"},
        },
        "scaling": {"method": "divide", "factor": 10},
        "output_suffix": "",
    },
    "AST_09T": {
        "name": "TIR Calibrated Radiance",
        "science_bands": {
            "SurfaceRadiance_B10": {"aster_band": 10, "group": "SurfaceRadiance", "description": "Surface Radiance Band 10"},
            "SurfaceRadiance_B11": {"aster_band": 11, "group": "SurfaceRadiance", "description": "Surface Radiance Band 11"},
            "SurfaceRadiance_B12": {"aster_band": 12, "group": "SurfaceRadiance", "description": "Surface Radiance Band 12"},
            "SurfaceRadiance_B13": {"aster_band": 13, "group": "SurfaceRadiance", "description": "Surface Radiance Band 13"},
            "SurfaceRadiance_B14": {"aster_band": 14, "group": "SurfaceRadiance", "description": "Surface Radiance Band 14"},
            "SkyIrradiance_B10":   {"aster_band": 10, "group": "SkyIrradiance",   "description": "Sky Irradiance Band 10"},
            "SkyIrradiance_B11":   {"aster_band": 11, "group": "SkyIrradiance",   "description": "Sky Irradiance Band 11"},
            "SkyIrradiance_B12":   {"aster_band": 12, "group": "SkyIrradiance",   "description": "Sky Irradiance Band 12"},
            "SkyIrradiance_B13":   {"aster_band": 13, "group": "SkyIrradiance",   "description": "Sky Irradiance Band 13"},
            "SkyIrradiance_B14":   {"aster_band": 14, "group": "SkyIrradiance",   "description": "Sky Irradiance Band 14"},
        },
        "scaling": {"method": "per_band_multiply", "factors": TIR_UCC},
        "output_suffix": "",   # SRA / SIR suffixes added per-group below
    },
    "AST_L1T": {
        "name": "Precision Terrain Corrected Radiance",
        "science_bands": {
            "VNIR_B01": {"aster_band": 1,  "group": "VNIR", "description": "VNIR Band 1"},
            "VNIR_B02": {"aster_band": 2,  "group": "VNIR", "description": "VNIR Band 2"},
            "VNIR_B3N": {"aster_band": 3,  "group": "VNIR", "description": "VNIR Band 3N"},
            "SWIR_B04": {"aster_band": 4,  "group": "SWIR", "description": "SWIR Band 4",  "optional": True},
            "SWIR_B05": {"aster_band": 5,  "group": "SWIR", "description": "SWIR Band 5",  "optional": True},
            "SWIR_B06": {"aster_band": 6,  "group": "SWIR", "description": "SWIR Band 6",  "optional": True},
            "SWIR_B07": {"aster_band": 7,  "group": "SWIR", "description": "SWIR Band 7",  "optional": True},
            "SWIR_B08": {"aster_band": 8,  "group": "SWIR", "description": "SWIR Band 8",  "optional": True},
            "SWIR_B09": {"aster_band": 9,  "group": "SWIR", "description": "SWIR Band 9",  "optional": True},
            "TIR_B10":  {"aster_band": 10, "group": "TIR",  "description": "TIR Band 10"},
            "TIR_B11":  {"aster_band": 11, "group": "TIR",  "description": "TIR Band 11"},
            "TIR_B12":  {"aster_band": 12, "group": "TIR",  "description": "TIR Band 12"},
            "TIR_B13":  {"aster_band": 13, "group": "TIR",  "description": "TIR Band 13"},
            "TIR_B14":  {"aster_band": 14, "group": "TIR",  "description": "TIR Band 14"},
        },
        # VNIR/SWIR: no scaling (raw DN).  TIR: (DN-1) × UCC handled specially in _process_granule.
        "scaling": None,
        "output_suffix": "",   # VNIR / SWIR / TIR suffixes added per-group below
    },
}

# Nodata sentinel values used in ASTER V004 GeoTIFFs
_NODATA_VALUES = {np.float32(v) for v in (-9999, 0, -32768, 32767, 65535, -1)}


# =============================================================================
# UTILITY FUNCTIONS  (ported from ASTER Data Toolkit v1.2)
# =============================================================================

def _is_qa_file(filepath):
    fu = os.path.basename(filepath).upper()
    if "QA_DATAPLANE" in fu or "QA_DATA_PLANE" in fu:
        return True
    if os.path.splitext(fu)[0].endswith("_QA"):
        return True
    return False


def _get_raster_files(folder, exclude_qa=False):
    files = set()
    for pat in ("*.tif", "*.tiff", "*.TIF", "*.TIFF", "*.dat", "*.DAT"):
        for f in glob.glob(os.path.join(folder, pat)):
            if exclude_qa and _is_qa_file(f):
                continue
            files.add(os.path.normpath(f))
    return sorted(files)


def _extract_granule_id(filepath):
    """Return the granule ID prefix from an ASTER V004 filename."""
    fname = os.path.splitext(os.path.basename(filepath))[0]
    m = re.match(r"(AST_(?:L1T|\w+?)_\d+_\d+)_", fname)
    return m.group(1) if m else None


def _group_files_by_granule(file_list, exclude_qa=False):
    groups = {}
    for fp in file_list:
        if exclude_qa and _is_qa_file(fp):
            continue
        gid = _extract_granule_id(fp)
        if gid:
            groups.setdefault(gid, []).append(fp)
    return groups


def _discover_granule_folders(parent):
    """Return list of (subfolder_path, folder_name, tif_list) for each granule subfolder."""
    result = []
    for entry in sorted(os.listdir(parent)):
        sf = os.path.join(parent, entry)
        if os.path.isdir(sf):
            tifs = _get_raster_files(sf)
            if tifs:
                result.append((sf, entry, tifs))
    return result


def _safe_output_path(path, overwrite=True):
    """
    If overwrite is True, return path as-is (will overwrite existing).
    If overwrite is False, auto-increment (_v2, _v3, …) to avoid overwriting.
    """
    if overwrite or not os.path.exists(path):
        return path
    base, ext = os.path.splitext(path)
    v = 2
    while True:
        candidate = f"{base}_v{v}{ext}"
        if not os.path.exists(candidate):
            return candidate
        v += 1


def _mask_nodata(data, nodata_val=None):
    """Convert nodata sentinels to NaN in a float32 array."""
    data = data.astype(np.float32)
    if nodata_val is not None:
        data[data == np.float32(nodata_val)] = np.nan
    for fv in _NODATA_VALUES:
        data[data == fv] = np.nan
    return data


def _classify_files(file_list, product_key):
    """
    Sort a list of single-band GeoTIFFs into science bands based on filename patterns.

    Returns:
        science_files : list of (filepath, band_key, band_info_dict)
        qa_files      : list of (filepath, label_str)
    """
    science, qa = [], []

    for fpath in file_list:
        fu = os.path.basename(fpath).upper()

        # ---- QA files ----
        if "QA_DATAPLANE" in fu or "QA_DATA_PLANE" in fu:
            if product_key == "AST_09T":
                label = ("SkyIrradiance_QA" if "_SIR_" in fu
                         else "SurfaceRadiance_QA" if "_SRA_" in fu
                         else "QA")
            else:
                label = "QA"
            qa.append((fpath, label))
            continue

        classified = False

        if product_key == "AST_05":
            m = re.search(r"_SRE_TIR_B(\d+)", fu)
            if m:
                bk = f"Band_{int(m.group(1))}"
                if bk in ASTER_PRODUCTS["AST_05"]["science_bands"]:
                    science.append((fpath, bk, ASTER_PRODUCTS["AST_05"]["science_bands"][bk]))
                    classified = True

        elif product_key == "AST_08":
            if "_SKT" in fu and "QA" not in fu:
                science.append((fpath, "Temperature",
                                ASTER_PRODUCTS["AST_08"]["science_bands"]["Temperature"]))
                classified = True

        elif product_key == "AST_09T":
            for prefix, grp in (("_SRA_TIR_B", "SurfaceRadiance"),
                                 ("_SIR_TIR_B", "SkyIrradiance")):
                m = re.search(prefix + r"(\d+)", fu)
                if m:
                    bk = f"{grp}_B{int(m.group(1))}"
                    if bk in ASTER_PRODUCTS["AST_09T"]["science_bands"]:
                        science.append((fpath, bk,
                                        ASTER_PRODUCTS["AST_09T"]["science_bands"][bk]))
                        classified = True

        elif product_key == "AST_07X":
            m_v = re.search(r"_SRF_VNIR_B(\d+N?)", fu)
            m_s = re.search(r"_SRF_SWIR_B(\d+)", fu)
            if m_v:
                bid = m_v.group(1)
                bk = ("Band_3N" if bid in ("3N", "03N")
                      else f"Band_{int(bid):02d}")
                if bk in ASTER_PRODUCTS["AST_07X"]["science_bands"]:
                    science.append((fpath, bk,
                                    ASTER_PRODUCTS["AST_07X"]["science_bands"][bk]))
                    classified = True
            elif m_s:
                bk = f"Band_{int(m_s.group(1)):02d}"
                if bk in ASTER_PRODUCTS["AST_07X"]["science_bands"]:
                    science.append((fpath, bk,
                                    ASTER_PRODUCTS["AST_07X"]["science_bands"][bk]))
                    classified = True

        elif product_key == "AST_L1T":
            m_vnir = re.search(r"_VNIR_B(\d+N?)", fu)
            m_swir = re.search(r"_SWIR_B(\d+)", fu)
            m_tir  = re.search(r"_TIR_B(\d+)", fu)
            if m_vnir:
                bid = m_vnir.group(1)
                bk = ("VNIR_B3N" if bid in ("3N", "03N")
                      else f"VNIR_B{int(bid):02d}")
                if bk in ASTER_PRODUCTS["AST_L1T"]["science_bands"]:
                    science.append((fpath, bk,
                                    ASTER_PRODUCTS["AST_L1T"]["science_bands"][bk]))
                    classified = True
            elif m_swir:
                bk = f"SWIR_B{int(m_swir.group(1)):02d}"
                if bk in ASTER_PRODUCTS["AST_L1T"]["science_bands"]:
                    science.append((fpath, bk,
                                    ASTER_PRODUCTS["AST_L1T"]["science_bands"][bk]))
                    classified = True
            elif m_tir:
                bk = f"TIR_B{int(m_tir.group(1))}"
                if bk in ASTER_PRODUCTS["AST_L1T"]["science_bands"]:
                    science.append((fpath, bk,
                                    ASTER_PRODUCTS["AST_L1T"]["science_bands"][bk]))
                    classified = True

    # Sort science bands by (group order, band number)
    _group_order = {"SurfaceRadiance": 0, "SkyIrradiance": 1, "": 0}
    def _sort_key(item):
        g = item[2].get("group", "")
        ab = item[2].get("aster_band") or 0
        return (_group_order.get(g, 9), ab)

    science.sort(key=_sort_key)
    qa.sort(key=lambda x: x[0])
    return science, qa


# =============================================================================
# I/O HELPERS
# =============================================================================

def _apply_scaling(data, product_key, aster_band_num):
    """Apply product scaling in-place and return the array.
    If the product defines a valid_range, pixels outside that range are set to NaN.
    """
    prod = ASTER_PRODUCTS.get(product_key)
    if not prod or not prod.get("scaling"):
        return data
    sc = prod["scaling"]
    valid = ~np.isnan(data)
    if sc["method"] == "divide":
        data[valid] = data[valid] / sc["factor"]
    elif sc["method"] == "per_band_multiply" and aster_band_num is not None:
        ucc = sc["factors"].get(aster_band_num, 1.0)
        # AST_09T raw DNs: Radiance = (DN - 1) × UCC
        data[valid] = (data[valid] - 1.0) * ucc

    # Mask pixels outside valid range (e.g. reflectance/emissivity must be 0–1)
    vr = prod.get("valid_range")
    if vr is not None:
        lo, hi = vr
        out_of_range = ~np.isnan(data) & ((data < lo) | (data > hi))
        n_masked = np.count_nonzero(out_of_range)
        if n_masked > 0:
            data[out_of_range] = np.nan

    return data


def _write_geotiff(arrays, output_path, ref_meta, descriptions,
                   compress='deflate', overwrite=True):
    """Write a multiband float32 GeoTIFF."""
    output_path = _safe_output_path(output_path, overwrite=overwrite)
    meta = ref_meta.copy()
    meta.update({
        'driver':    'GTiff',
        'count':     len(arrays),
        'dtype':     'float32',
        'nodata':    np.nan,
        'compress':  compress,
        'predictor': 2,
        'zlevel':    6,
    })
    with rasterio.open(output_path, 'w', **meta) as dst:
        for i, (arr, desc) in enumerate(zip(arrays, descriptions), 1):
            dst.write(arr.astype(np.float32), i)
            dst.set_band_description(i, desc)
    return output_path


def _write_envi(arrays, output_base, descriptions, wavelengths, crs_wkt,
                transform, overwrite=True):
    """Write a multiband float32 ENVI .dat + .hdr pair."""
    dat_path = _safe_output_path(output_base + '.dat', overwrite=overwrite)
    hdr_path = os.path.splitext(dat_path)[0] + '.hdr'

    # ENVI nodata — use -9999 as the fill value (avoids confusion with valid 0 reflectance)
    envi_nodata = -9999.0

    rows, cols = arrays[0].shape
    meta = {
        'driver': 'ENVI',
        'height': rows,
        'width':  cols,
        'count':  len(arrays),
        'dtype':  'float32',
        'nodata': envi_nodata,
    }
    if crs_wkt:
        meta['crs'] = CRS.from_wkt(crs_wkt)
    if transform:
        meta['transform'] = transform

    with rasterio.open(dat_path, 'w', **meta) as dst:
        for i, (arr, desc) in enumerate(zip(arrays, descriptions), 1):
            out = arr.astype(np.float32).copy()
            out[np.isnan(out)] = envi_nodata
            dst.write(out, i)
            dst.set_band_description(i, desc)

    # Append wavelength / band-name metadata to the header
    # Also fix nodata if rasterio wrote 'nan' instead of our numeric value
    if os.path.exists(hdr_path):
        with open(hdr_path, 'r') as f:
            hdr = f.read()

        # Fix nodata: replace 'nan' with '-9999'
        hdr = re.sub(r'data ignore value\s*=\s*nan',
                      f'data ignore value = {envi_nodata}', hdr,
                      flags=re.IGNORECASE)

        extra = ''
        if 'wavelength =' not in hdr and wavelengths:
            extra += f"\nwavelength = {{{', '.join(f'{w:.4f}' for w in wavelengths)}}}\n"
        if 'wavelength units =' not in hdr:
            extra += 'wavelength units = Micrometers\n'
        if 'band names =' not in hdr and descriptions:
            extra += f"band names = {{{', '.join(descriptions)}}}\n"

        # Rewrite the entire header (to apply the nodata fix + extras)
        with open(hdr_path, 'w') as f:
            f.write(hdr)
            if extra:
                f.write(extra)

    # Remove spurious .aux.xml that rasterio's ENVI driver creates —
    # it contains no geolocation and can confuse ENVI
    for aux_ext in ('.dat.aux.xml', '.aux.xml'):
        spurious = os.path.splitext(dat_path)[0] + aux_ext
        if os.path.isfile(spurious):
            os.remove(spurious)
    # Also check for dat_path + .aux.xml
    if os.path.isfile(dat_path + '.aux.xml'):
        os.remove(dat_path + '.aux.xml')

    return dat_path


# =============================================================================
# GRANULE PROCESSOR
# =============================================================================

def _process_granule(granule_key, science_files, product_key,
                     out_tif_dir, out_envi_dir,
                     write_geotiff, write_envi, overwrite=True,
                     qa_files=None):
    """
    Stack and scale all science bands for one granule, then write outputs.
    For AST_09T, writes two separate files (SurfaceRadiance, SkyIrradiance).
    If qa_files are provided, writes a separate QA multiband file.
    """
    if not science_files:
        print(f"  [skip] {granule_key} — no science files classified")
        return 0

    # Determine output groups
    if product_key == "AST_09T":
        groups = [
            ('SRA', [x for x in science_files if x[2].get('group') == 'SurfaceRadiance']),
            ('SIR', [x for x in science_files if x[2].get('group') == 'SkyIrradiance']),
        ]
    elif product_key == "AST_L1T":
        groups = [
            ('VNIR', [x for x in science_files if x[2].get('group') == 'VNIR']),
            ('SWIR', [x for x in science_files if x[2].get('group') == 'SWIR']),
            ('TIR',  [x for x in science_files if x[2].get('group') == 'TIR']),
        ]
    else:
        groups = [('', science_files)]

    files_written = 0

    for suffix, band_list in groups:
        if not band_list:
            continue

        # Read + mask + scale each band
        arrays, descriptions, wavelengths = [], [], []
        ref_meta, crs_wkt, transform = None, None, None

        # Is this the L1T TIR group? Needs special DN → radiance conversion.
        is_l1t_tir = (product_key == "AST_L1T" and suffix == "TIR")

        for fpath, band_key, band_info in band_list:
            with rasterio.open(fpath) as src:
                data = _mask_nodata(src.read(1), src.nodata)
                if ref_meta is None:
                    ref_meta = src.meta.copy()
                if crs_wkt is None and src.crs:
                    crs_wkt = src.crs.to_wkt()
                if transform is None:
                    transform = src.transform

            abn = band_info.get('aster_band')

            if is_l1t_tir:
                # Radiance = (DN - 1) × UCC
                valid = ~np.isnan(data)
                ucc = TIR_UCC.get(abn, 1.0)
                data[valid] = (data[valid] - 1.0) * ucc
            else:
                # All other products/groups: use standard scaling (or no-op if scaling=None)
                data = _apply_scaling(data, product_key, abn)

            arrays.append(data)
            descriptions.append(band_info['description'])
            wavelengths.append(ASTER_WAVELENGTHS.get(abn, 0.0) if abn else 0.0)

        # Build output base name
        gid = _extract_granule_id(band_list[0][0]) or granule_key
        base_name = f"{gid}_{suffix}" if suffix else gid

        # Write GeoTIFF
        if write_geotiff and out_tif_dir:
            out_path = os.path.join(out_tif_dir, f"{base_name}.tif")
            result = _write_geotiff(arrays, out_path, ref_meta, descriptions,
                                    overwrite=overwrite)
            print(f"  [GeoTIFF] {os.path.basename(result)}")
            files_written += 1

        # Write ENVI
        if write_envi and out_envi_dir:
            out_base = os.path.join(out_envi_dir, base_name)
            result = _write_envi(arrays, out_base, descriptions,
                                 wavelengths, crs_wkt, transform,
                                 overwrite=overwrite)
            print(f"  [ENVI]    {os.path.basename(result)}")
            files_written += 1

    # --- Write QA file(s) ---
    if qa_files:
        # Group QA files by label for AST_09T (SkyIrradiance_QA, SurfaceRadiance_QA)
        qa_groups = {}
        for fpath, label in qa_files:
            qa_groups.setdefault(label, []).append(fpath)

        for qa_label, qa_paths in sorted(qa_groups.items()):
            qa_arrays, qa_descriptions = [], []
            qa_ref_meta, qa_crs_wkt, qa_transform = None, None, None

            for fpath in sorted(qa_paths):
                with rasterio.open(fpath) as src:
                    data = src.read(1).astype(np.float32)
                    if qa_ref_meta is None:
                        qa_ref_meta = src.meta.copy()
                    if qa_crs_wkt is None and src.crs:
                        qa_crs_wkt = src.crs.to_wkt()
                    if qa_transform is None:
                        qa_transform = src.transform

                qa_arrays.append(data)
                # Build description from filename
                fname = os.path.basename(fpath)
                if 'DATAPLANE2' in fname.upper():
                    qa_descriptions.append(f'{qa_label} DataPlane2')
                else:
                    qa_descriptions.append(f'{qa_label} DataPlane')

            if not qa_arrays or qa_ref_meta is None:
                continue

            gid = _extract_granule_id(qa_paths[0]) or granule_key
            if qa_label == 'QA':
                qa_base_name = f"{gid}_QA"
            else:
                qa_base_name = f"{gid}_{qa_label}"

            if write_geotiff and out_tif_dir:
                out_path = os.path.join(out_tif_dir, f"{qa_base_name}.tif")
                result = _write_geotiff(qa_arrays, out_path, qa_ref_meta,
                                        qa_descriptions, overwrite=overwrite)
                print(f"  [GeoTIFF] {os.path.basename(result)}  (QA)")
                files_written += 1

            if write_envi and out_envi_dir:
                out_base = os.path.join(out_envi_dir, qa_base_name)
                result = _write_envi(qa_arrays, out_base, qa_descriptions,
                                     [0.0] * len(qa_arrays), qa_crs_wkt,
                                     qa_transform, overwrite=overwrite)
                print(f"  [ENVI]    {os.path.basename(result)}  (QA)")
                files_written += 1

    return files_written


# =============================================================================
# PUBLIC ENTRY POINT
# =============================================================================

def process_directory(input_dir, output_dir=None,
                      aster_product='AST_05',
                      write_geotiff=True, write_envi=True,
                      overwrite=True):
    """
    Convert all ASTER granules in input_dir to analysis-ready multiband files.

    Parameters
    ----------
    input_dir    : str   Path to folder containing granule subfolders or flat TIFFs.
    output_dir   : str   Output root.  Defaults to <input_dir>/converted/ASTER.
    aster_product: str   One of: 'AST_05', 'AST_07', 'AST_08', 'AST_09T', 'AST_L1T'.
                         'AST_07' and 'AST_07XT' are treated identically.
    write_geotiff: bool  Write GeoTIFF output.
    write_envi   : bool  Write ENVI output.
    overwrite    : bool  If True, overwrite existing files. If False, auto-increment.
    """

    # Normalise product key
    if aster_product in ('AST_07', 'AST_07XT'):
        product_key = 'AST_07X'
        display_name = 'AST_07 / AST_07XT'
    else:
        product_key = aster_product
        display_name = aster_product

    if product_key not in ASTER_PRODUCTS:
        print(f"ERROR: Unknown ASTER product '{aster_product}'.  "
              f"Choose from: AST_05, AST_07, AST_08, AST_09T, AST_L1T")
        return

    if not write_geotiff and not write_envi:
        print("ERROR: At least one of write_geotiff or write_envi must be True.")
        return

    # Set up output directories
    if output_dir is None:
        output_dir = os.path.join(input_dir, 'converted', 'ASTER')

    out_tif_dir  = os.path.join(output_dir, 'GeoTIFF') if write_geotiff else None
    out_envi_dir = os.path.join(output_dir, 'ENVI')    if write_envi    else None

    for d in filter(None, (out_tif_dir, out_envi_dir)):
        os.makedirs(d, exist_ok=True)

    print(f"ASTER Converter  |  product: {display_name}")
    print(f"Input  : {input_dir}")
    print(f"Output : {output_dir}")
    print(f"Format : {'GeoTIFF ' if write_geotiff else ''}{'ENVI' if write_envi else ''}")
    print("-" * 60)

    # ------------------------------------------------------------------
    # Discover granules — subfolder mode first, then flat-folder fallback
    # ------------------------------------------------------------------
    granule_folders = _discover_granule_folders(input_dir)

    if granule_folders:
        print(f"Found {len(granule_folders)} granule subfolder(s)  [subfolder mode]")
        granules = {}
        for sf, folder_name, tifs in granule_folders:
            sci, qa = _classify_files(tifs, product_key)
            if sci:
                granules[folder_name] = (sci, qa)
            else:
                print(f"  [skip] {folder_name} — no {display_name} files recognised")
    else:
        print("No subfolders found — scanning flat folder...")
        all_tifs = _get_raster_files(input_dir)
        if not all_tifs:
            print("ERROR: No GeoTIFF/ENVI files found in input directory.")
            return
        groups = _group_files_by_granule(all_tifs)
        if not groups:
            print("ERROR: Could not identify any ASTER granule IDs from filenames.")
            return
        print(f"Found {len(groups)} granule(s)  [flat-folder mode]")
        granules = {}
        for gid, files in sorted(groups.items()):
            sci, qa = _classify_files(files, product_key)
            if sci:
                granules[gid] = (sci, qa)
            else:
                print(f"  [skip] {gid} — no {display_name} files recognised")

    if not granules:
        print(f"No {display_name} files found.  "
              f"Ensure the correct product is selected and files are present.")
        return

    print(f"\nProcessing {len(granules)} granule(s)...")
    print("-" * 60)

    total_files = 0
    for i, (gkey, (sci, qa)) in enumerate(granules.items(), 1):
        print(f"[{i}/{len(granules)}] {gkey}  ({len(sci)} band(s), {len(qa)} QA file(s))")
        n = _process_granule(gkey, sci, product_key,
                             out_tif_dir, out_envi_dir,
                             write_geotiff, write_envi,
                             overwrite=overwrite,
                             qa_files=qa if qa else None)
        total_files += n

    print("-" * 60)
    print(f"Done.  {total_files} file(s) written to: {output_dir}")


# =============================================================================
# V003 HDF PRODUCT DEFINITIONS
# =============================================================================

V003_PRODUCTS = {
    "AST_05": {
        "name": "TIR Surface Emissivity (V003 HDF)",
        "bands": {
            "Band10": {"aster_band": 10, "description": "Surface Emissivity Band 10"},
            "Band11": {"aster_band": 11, "description": "Surface Emissivity Band 11"},
            "Band12": {"aster_band": 12, "description": "Surface Emissivity Band 12"},
            "Band13": {"aster_band": 13, "description": "Surface Emissivity Band 13"},
            "Band14": {"aster_band": 14, "description": "Surface Emissivity Band 14"},
        },
        "scaling": {"method": "multiply", "factor": 0.001},
        "valid_range": (0.0, 1.0),
        "dtype": "int32",
        "shape": (700, 830),
        "has_aux": True,
    },
    "AST_07": {
        "name": "Surface Reflectance VNIR/SWIR (V003 HDF)",
        "bands": {
            "Band1":  {"aster_band": 1, "description": "VNIR Reflectance Band 1"},
            "Band2":  {"aster_band": 2, "description": "VNIR Reflectance Band 2"},
            "Band3N": {"aster_band": 3, "description": "VNIR Reflectance Band 3N"},
            "Band4":  {"aster_band": 4, "description": "SWIR Reflectance Band 4",  "optional": True},
            "Band5":  {"aster_band": 5, "description": "SWIR Reflectance Band 5",  "optional": True},
            "Band6":  {"aster_band": 6, "description": "SWIR Reflectance Band 6",  "optional": True},
            "Band7":  {"aster_band": 7, "description": "SWIR Reflectance Band 7",  "optional": True},
            "Band8":  {"aster_band": 8, "description": "SWIR Reflectance Band 8",  "optional": True},
            "Band9":  {"aster_band": 9, "description": "SWIR Reflectance Band 9",  "optional": True},
        },
        "scaling": {"method": "multiply", "factor": 0.001},
        "valid_range": (0.0, 1.0),
        "dtype": "int16",
        "shape": None,   # varies: 4200×4980 (VNIR) or 2100×2490 (SWIR)
        "has_aux": False,
    },
    "AST_08": {
        "name": "Surface Kinetic Temperature (V003 HDF)",
        "bands": {
            "KineticTemperature": {"aster_band": None, "description": "Surface Kinetic Temperature"},
        },
        "scaling": {"method": "multiply", "factor": 0.1},
        "dtype": "int32",
        "shape": (700, 830),
        "has_aux": False,
    },
    "AST_09T": {
        "name": "TIR Calibrated Radiance (V003 HDF)",
        "bands": {
            "Band10": {"aster_band": 10, "description": "Surface Radiance Band 10"},
            "Band11": {"aster_band": 11, "description": "Surface Radiance Band 11"},
            "Band12": {"aster_band": 12, "description": "Surface Radiance Band 12"},
            "Band13": {"aster_band": 13, "description": "Surface Radiance Band 13"},
            "Band14": {"aster_band": 14, "description": "Surface Radiance Band 14"},
        },
        "scaling": {"method": "per_band_multiply", "factors": TIR_UCC},
        "dtype": "int16",
        "shape": (700, 830),
        "has_aux": True,
    },
}

# Nodata values used in V003 HDF products
_V003_NODATA = {-9999, 0, -32768, 32767}


def _extract_v003_acquisition_id(hdf_path):
    """
    Extract the acquisition timestamp portion from a V003 filename to
    match VNIR/SWIR pairs.

    Filename patterns:
        AST_07_003MMDDYYYYHHMMSS_YYYYMMDDHHMMSS_NNNNN.hdf
        AST_07XT_003MMDDYYYYHHMMSS_YYYYMMDDHHMMSS_NNNNN.hdf

    Returns the '003MMDDYYYYHHMMSS' acquisition portion only (not processing time).
    """
    fname = os.path.splitext(os.path.basename(hdf_path))[0]
    # Strip the product prefix to get the acquisition timestamp
    # Pattern: AST_07XT_003MMDDYYYYHHMMSS_YYYYMMDDHHMMSS_NNNNN
    #                    ^^^^^^^^^^^^^^^^^ this part (observation time)
    m = re.match(r'AST_07(?:XT)?_(003\d{14})', fname, re.IGNORECASE)
    if m:
        return m.group(1)
    return None


def _detect_v003_07_type(hdf_path):
    """
    Determine if a V003 AST_07/07XT file contains VNIR or SWIR bands
    by reading its SDS names.  Returns 'vnir', 'swir', or None.
    """
    if not _HAS_PYHDF:
        return None
    try:
        sd = SD(hdf_path, SDC.READ)
        sds_names = list(sd.datasets().keys())
        sd.end()

        has_vnir = any(n in sds_names for n in ('Band1', 'Band2', 'Band3N'))
        has_swir = any(n in sds_names for n in ('Band4', 'Band5', 'Band6',
                                                 'Band7', 'Band8', 'Band9'))
        if has_vnir and not has_swir:
            return 'vnir'
        elif has_swir and not has_vnir:
            return 'swir'
        elif has_vnir and has_swir:
            return 'both'
        return None
    except Exception:
        return None


def _downsample_2x2(arr):
    """Downsample a 2D array by 2× using block averaging, preserving NaN."""
    rows, cols = arr.shape
    # Trim to even dimensions
    rows_even = rows - (rows % 2)
    cols_even = cols - (cols % 2)
    trimmed = arr[:rows_even, :cols_even]

    # Reshape into 2×2 blocks
    blocks = trimmed.reshape(rows_even // 2, 2, cols_even // 2, 2)

    # nanmean to preserve radiometry while handling NaN fill values
    with np.errstate(all='ignore'):
        result = np.nanmean(blocks, axis=(1, 3))

    return result


# =============================================================================
# V003 GEOLOCATION
# =============================================================================

def _parse_met_file(met_path):
    """
    Parse an ASTER V003 .met sidecar file and extract corner coordinates
    and UTM zone.  Supports two formats:

    Format A (embedded in .aux.xml or some .met files):
        UPPERLEFT = "lat, lon"
        UTMZONECODE = 18

    Format B (NASA ODL .met files):
        GRINGPOINTLATITUDE = (lat1, lat2, lat3, lat4)
        GRINGPOINTLONGITUDE = (lon1, lon2, lon3, lon4)
        GRINGPOINTSEQUENCENO = (1, 2, 3, 4)
        Sequence: 1=UL, 2=UR, 3=LR, 4=LL

    Also checks for ASTERMapProjection / UTM info in the .met.

    Returns dict with keys:
        'upper_left', 'upper_right', 'lower_left', 'lower_right'
            — each a (lat, lon) tuple
        'utm_zone' — integer UTM zone code (or None)
    """
    corners = {}
    utm_zone = None

    with open(met_path, 'r') as f:
        content = f.read()

    # --- Format A: UPPERLEFT/LOWERRIGHT style ---
    for tag, key in (('UPPERLEFT',  'upper_left'),
                     ('UPPERRIGHT', 'upper_right'),
                     ('LOWERLEFT',  'lower_left'),
                     ('LOWERRIGHT', 'lower_right')):
        m = re.search(tag + r'\s*=\s*"?\s*([-\d.]+)\s*,\s*([-\d.]+)\s*"?', content)
        if m:
            corners[key] = (float(m.group(1)), float(m.group(2)))

    # UTM zone (Format A style)
    m = re.search(r'UTMZONECODE\s*=\s*(\d+)', content)
    if m:
        utm_zone = int(m.group(1))

    if len(corners) == 4:
        return corners, utm_zone

    # --- Format B: NASA ODL with GRING polygon ---
    # Extract lat tuple: VALUE = (-15.58, -15.68, -16.25, -16.14)
    m_lat = re.search(
        r'GRINGPOINTLATITUDE.*?VALUE\s*=\s*\(\s*([-\d.,\s]+)\s*\)',
        content, re.DOTALL)
    m_lon = re.search(
        r'GRINGPOINTLONGITUDE.*?VALUE\s*=\s*\(\s*([-\d.,\s]+)\s*\)',
        content, re.DOTALL)

    if m_lat and m_lon:
        try:
            lats = [float(v.strip()) for v in m_lat.group(1).split(',')]
            lons = [float(v.strip()) for v in m_lon.group(1).split(',')]

            if len(lats) >= 4 and len(lons) >= 4:
                # ASTER GRING convention: sequence 1=UL, 2=UR, 3=LR, 4=LL
                corners = {
                    'upper_left':  (lats[0], lons[0]),
                    'upper_right': (lats[1], lons[1]),
                    'lower_right': (lats[2], lons[2]),
                    'lower_left':  (lats[3], lons[3]),
                }
        except (ValueError, IndexError):
            pass

    # Check for UTM in ODL format
    if utm_zone is None:
        if 'Universal Transverse Mercator' in content or '"UTM"' in content:
            # Try to infer UTM zone from central meridian or scene center
            # The .met may not have explicit zone, but we can compute from lon
            if corners:
                avg_lon = sum(c[1] for c in corners.values()) / len(corners)
                utm_zone = int((avg_lon + 180) / 6) + 1

    return corners, utm_zone


def _parse_aux_xml(aux_path):
    """
    Parse an ASTER V003 .aux.xml (PAMDataset format) to extract corner
    coordinates, UTM zone, and CRS WKT.

    The .aux.xml contains <MDI key="UPPERLEFT">lat, lon</MDI> etc. in the
    per-subdataset metadata, plus a <Metadata domain="GEOLOCATION"> block
    with an SRS WKT string.

    Returns (corners_dict, utm_zone, crs_wkt) or ({}, None, None).
    """
    corners = {}
    utm_zone = None
    crs_wkt = None

    try:
        with open(aux_path, 'r') as f:
            content = f.read()

        # Extract corner coordinates from <MDI> keys
        for tag, key in (('UPPERLEFT',  'upper_left'),
                         ('UPPERRIGHT', 'upper_right'),
                         ('LOWERLEFT',  'lower_left'),
                         ('LOWERRIGHT', 'lower_right')):
            m = re.search(
                r'<MDI\s+key="' + tag + r'">\s*([-\d.]+)\s*,\s*([-\d.]+)\s*</MDI>',
                content)
            if m:
                corners[key] = (float(m.group(1)), float(m.group(2)))

        # Extract UTM zone from UTMZONECODE (any band-specific variant)
        m = re.search(r'<MDI\s+key="UTMZONECODE\d*">\s*(\d+)\s*</MDI>', content)
        if m:
            utm_zone = int(m.group(1))

        # Extract SRS WKT from GEOLOCATION metadata domain
        m = re.search(r'<MDI\s+key="SRS">(.*?)</MDI>', content, re.DOTALL)
        if m:
            crs_wkt = m.group(1).strip()

    except Exception:
        pass

    return corners, utm_zone, crs_wkt


def _parse_hdf_metadata(hdf_path):
    """
    Read corner coordinates and UTM zone directly from the HDF file's
    global/SDS attributes using pyhdf.  This works even without sidecar files.

    Returns (corners_dict, utm_zone) or ({}, None).
    """
    if not _HAS_PYHDF:
        return {}, None

    corners = {}
    utm_zone = None

    try:
        sd = SD(hdf_path, SDC.READ)
        attrs = sd.attributes()

        # Corner coordinates — stored as "lat, lon" strings
        for tag, key in (('UPPERLEFT',  'upper_left'),
                         ('UPPERRIGHT', 'upper_right'),
                         ('LOWERLEFT',  'lower_left'),
                         ('LOWERRIGHT', 'lower_right')):
            val = attrs.get(tag)
            if val and isinstance(val, str):
                parts = val.split(',')
                if len(parts) == 2:
                    try:
                        corners[key] = (float(parts[0].strip()),
                                        float(parts[1].strip()))
                    except ValueError:
                        pass

        # UTM zone — try various key patterns
        for key_name in ('UTMZONECODE', 'UTMZONECODE1', 'UTMZONECODE10'):
            val = attrs.get(key_name)
            if val is not None:
                try:
                    utm_zone = int(val)
                    break
                except (ValueError, TypeError):
                    pass

        sd.end()
    except Exception:
        pass

    return corners, utm_zone


def _build_geotransform_from_corners(corners, rows, cols, utm_zone=None):
    """
    Build a rasterio Affine transform and CRS from corner coordinates.

    Uses all four corners to compute a proper affine transform that accounts
    for rotation (ASTER swaths are typically rotated ~8-10° from north).

    Parameters
    ----------
    corners : dict with 'upper_left', 'upper_right', 'lower_left', 'lower_right'
              Each value is (lat, lon).
    rows, cols : int  Image dimensions.
    utm_zone : int or None.  If provided, builds a UTM CRS.

    Returns
    -------
    transform : rasterio.transform.Affine
    crs : rasterio.crs.CRS
    """
    ul_lat, ul_lon = corners['upper_left']
    ur_lat, ur_lon = corners['upper_right']
    ll_lat, ll_lon = corners['lower_left']
    lr_lat, lr_lon = corners['lower_right']

    if utm_zone is not None and utm_zone > 0:
        # Determine hemisphere from average latitude
        avg_lat = (ul_lat + ur_lat + ll_lat + lr_lat) / 4.0
        epsg = 32600 + utm_zone if avg_lat >= 0 else 32700 + utm_zone
        crs = CRS.from_epsg(epsg)

        try:
            from pyproj import Transformer
            transformer = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}",
                                                always_xy=True)
            ul_x, ul_y = transformer.transform(ul_lon, ul_lat)
            ur_x, ur_y = transformer.transform(ur_lon, ur_lat)
            ll_x, ll_y = transformer.transform(ll_lon, ll_lat)
            lr_x, lr_y = transformer.transform(lr_lon, lr_lat)
        except ImportError:
            print("  WARNING: pyproj not available — falling back to geographic CRS")
            crs = CRS.from_epsg(4326)
            ul_x, ul_y = ul_lon, ul_lat
            ur_x, ur_y = ur_lon, ur_lat
            ll_x, ll_y = ll_lon, ll_lat
            lr_x, lr_y = lr_lon, lr_lat
    else:
        crs = CRS.from_epsg(4326)
        ul_x, ul_y = ul_lon, ul_lat
        ur_x, ur_y = ur_lon, ur_lat
        ll_x, ll_y = ll_lon, ll_lat
        lr_x, lr_y = lr_lon, lr_lat

    # Build affine transform from corner points using least-squares fit.
    # The affine maps pixel (col, row) -> (x, y):
    #   x = a*col + b*row + c
    #   y = d*col + e*row + f
    # We have 4 corner points: (0,0)=UL, (cols,0)=UR, (0,rows)=LL, (cols,rows)=LR
    pixel_coords = np.array([
        [0,    0],       # UL
        [cols,  0],      # UR
        [0,    rows],    # LL
        [cols,  rows],   # LR
    ], dtype=np.float64)

    x_coords = np.array([ul_x, ur_x, ll_x, lr_x])
    y_coords = np.array([ul_y, ur_y, ll_y, lr_y])

    # Solve for affine coefficients: [col, row, 1] -> x  and  [col, row, 1] -> y
    A = np.column_stack([pixel_coords, np.ones(4)])
    # Least-squares solve
    x_coeffs, _, _, _ = np.linalg.lstsq(A, x_coords, rcond=None)
    y_coeffs, _, _, _ = np.linalg.lstsq(A, y_coords, rcond=None)

    # Affine(a, b, c, d, e, f) where:
    #   x = a*col + b*row + c
    #   y = d*col + e*row + f
    transform = rasterio.transform.Affine(
        x_coeffs[0], x_coeffs[1], x_coeffs[2],
        y_coeffs[0], y_coeffs[1], y_coeffs[2])

    return transform, crs


def _build_geotransform_from_tiepoints(hdf_path, rows, cols, utm_zone=None):
    """
    Build a geotransform by reading the Longitude/GeodeticLatitude tie-point
    arrays from the HDF file and computing an affine fit.

    The tie-point grids are 11×11 with LINE_STEP=70, PIXEL_STEP=83.
    We use the corner tie-points (which correspond to the actual image corners)
    to build the affine transform.

    Returns (transform, crs) or (None, None).
    """
    if not _HAS_PYHDF:
        return None, None

    try:
        sd = SD(hdf_path, SDC.READ)
        sds_names = list(sd.datasets().keys())

        if 'Longitude' not in sds_names or 'GeodeticLatitude' not in sds_names:
            sd.end()
            return None, None

        lon_data = sd.select('Longitude').get()
        lat_data = sd.select('GeodeticLatitude').get()
        sd.end()

        # Extract corner coordinates from the tie-point grid
        corners = {
            'upper_left':  (float(lat_data[0, 0]),   float(lon_data[0, 0])),
            'upper_right': (float(lat_data[0, -1]),   float(lon_data[0, -1])),
            'lower_left':  (float(lat_data[-1, 0]),   float(lon_data[-1, 0])),
            'lower_right': (float(lat_data[-1, -1]),  float(lon_data[-1, -1])),
        }

        # Determine UTM zone from the data if not provided
        if utm_zone is None:
            avg_lon = (lon_data[0, 0] + lon_data[0, -1] +
                       lon_data[-1, 0] + lon_data[-1, -1]) / 4.0
            utm_zone = int((avg_lon + 180) / 6) + 1

        return _build_geotransform_from_corners(corners, rows, cols, utm_zone)

    except Exception as e:
        print(f"  WARNING: Could not read tie-point arrays: {e}")
        return None, None


def _find_met_file(hdf_path):
    """Find the .met sidecar file for a V003 HDF file."""
    base = os.path.splitext(hdf_path)[0]
    for ext in ('.met', '.hdf.met'):
        candidate = base + ext
        if os.path.isfile(candidate):
            return candidate
    # Also check for the full filename + .met
    candidate = hdf_path + '.met'
    if os.path.isfile(candidate):
        return candidate
    return None


def _find_aux_xml(hdf_path):
    """Find the .aux.xml sidecar file for a V003 HDF file."""
    for candidate in (hdf_path + '.aux.xml',
                      os.path.splitext(hdf_path)[0] + '.aux.xml'):
        if os.path.isfile(candidate):
            return candidate
    return None


# =============================================================================
# V003 HDF I/O
# =============================================================================

def _detect_v003_product(hdf_path):
    """
    Detect the ASTER V003 product type from the HDF filename.
    Returns product key (e.g. 'AST_05') or None.
    """
    fname = os.path.basename(hdf_path).upper()
    if fname.startswith('AST_09T'):
        return 'AST_09T'
    elif fname.startswith('AST_05'):
        return 'AST_05'
    elif fname.startswith('AST_07'):
        return 'AST_07'
    elif fname.startswith('AST_08'):
        return 'AST_08'
    return None


def _read_v003_hdf(hdf_path, product_key):
    """
    Read science bands from a V003 HDF file using pyhdf.

    Returns
    -------
    arrays : list of 2D numpy float32 arrays
    descriptions : list of str
    wavelengths : list of float
    shape : (rows, cols) of the first band
    """
    if not _HAS_PYHDF:
        raise ImportError(
            "pyhdf is required for V003 HDF files.  "
            "Install with: conda install -c conda-forge pyhdf")

    prod = V003_PRODUCTS[product_key]
    sd = SD(hdf_path, SDC.READ)

    # List all SDS names in the file
    sds_names = list(sd.datasets().keys())

    arrays = []
    descriptions = []
    wavelengths = []
    shape = None

    for sds_key, band_info in prod['bands'].items():
        # Find the matching SDS — exact name first, then case-insensitive
        matched_name = None
        if sds_key in sds_names:
            matched_name = sds_key
        else:
            # Try case-insensitive match
            sds_key_lower = sds_key.lower()
            for name in sds_names:
                if name.lower() == sds_key_lower:
                    matched_name = name
                    break

        if matched_name is None:
            if band_info.get('optional'):
                continue   # silently skip optional bands not present
            print(f"  WARNING: SDS '{sds_key}' not found in {os.path.basename(hdf_path)}")
            print(f"           Available SDS: {sds_names}")
            continue

        sds = sd.select(matched_name)
        data = sds.get().astype(np.float32)
        sds.endaccess()

        if shape is None:
            shape = data.shape

        # Mask nodata
        for nd in _V003_NODATA:
            data[data == float(nd)] = np.nan

        # Apply scaling
        valid = ~np.isnan(data)
        sc = prod['scaling']
        abn = band_info.get('aster_band')

        if sc['method'] == 'multiply':
            data[valid] = data[valid] * sc['factor']
        elif sc['method'] == 'per_band_multiply' and abn is not None:
            ucc = sc['factors'].get(abn, 1.0)
            data[valid] = data[valid] * ucc

        # Mask pixels outside valid range (e.g. reflectance/emissivity 0–1)
        vr = prod.get('valid_range')
        if vr is not None:
            lo, hi = vr
            out_of_range = ~np.isnan(data) & ((data < lo) | (data > hi))
            if np.any(out_of_range):
                data[out_of_range] = np.nan

        arrays.append(data)
        descriptions.append(band_info['description'])
        wavelengths.append(ASTER_WAVELENGTHS.get(abn, 0.0) if abn else 0.0)

    sd.end()
    return arrays, descriptions, wavelengths, shape


def _get_v003_geolocation(hdf_path, product_key, shape):
    """
    Determine the geolocation (transform + CRS) for a V003 HDF file.

    Strategy (in priority order):
        1. HDF tie-point arrays (Longitude/GeodeticLatitude) — most accurate
        2. .aux.xml — corner coords, UTM zone, SRS WKT
        3. .met file — corner coords (UPPERLEFT style or GRING polygon)
        4. HDF global attributes — UPPERLEFT/LOWERRIGHT keys

    Returns (transform, crs) or (None, None).
    """
    rows, cols = shape

    # Determine UTM zone from aux.xml or met (needed for tie-point approach)
    utm_zone = None
    crs_wkt = None

    prod = V003_PRODUCTS[product_key]
    if prod.get('has_aux'):
        aux_path = _find_aux_xml(hdf_path)
        if aux_path:
            aux_corners, utm_zone, crs_wkt = _parse_aux_xml(aux_path)
            if utm_zone:
                print(f"  UTM zone {utm_zone} from .aux.xml")

    if utm_zone is None:
        met_path = _find_met_file(hdf_path)
        if met_path:
            met_corners, utm_zone = _parse_met_file(met_path)
            if utm_zone:
                print(f"  UTM zone {utm_zone} from .met")

    # Strategy 1: Tie-point arrays from HDF (best — uses actual pixel coords)
    transform, crs = _build_geotransform_from_tiepoints(
        hdf_path, rows, cols, utm_zone)
    if transform is not None:
        # NOTE: We do NOT override CRS with .aux.xml WKT here because
        # the NASA metadata often says "Northern Hemisphere" even for
        # southern-hemisphere scenes.  The EPSG computed from the actual
        # latitude in _build_geotransform_from_corners is authoritative.
        print(f"  Geolocation: from HDF tie-point arrays ({crs})")
        return transform, crs

    # Strategy 2: .aux.xml corner coordinates
    if prod.get('has_aux'):
        aux_path = _find_aux_xml(hdf_path)
        if aux_path:
            corners, uz, wkt = _parse_aux_xml(aux_path)
            if len(corners) == 4:
                transform, crs = _build_geotransform_from_corners(
                    corners, rows, cols, uz)
                if transform is not None:
                    print(f"  Geolocation: from .aux.xml corners ({crs})")
                    return transform, crs

    # Strategy 3: .met file corner coordinates
    met_path = _find_met_file(hdf_path)
    if met_path:
        corners, uz = _parse_met_file(met_path)
        if len(corners) == 4:
            transform, crs = _build_geotransform_from_corners(
                corners, rows, cols, uz)
            if transform is not None:
                print(f"  Geolocation: from .met corners ({crs})")
                return transform, crs

    # Strategy 4: HDF global attributes
    corners, uz = _parse_hdf_metadata(hdf_path)
    if len(corners) == 4:
        transform, crs = _build_geotransform_from_corners(
            corners, rows, cols, uz)
        if transform is not None:
            print(f"  Geolocation: from HDF metadata ({crs})")
            return transform, crs

    print(f"  WARNING: No geolocation found for {os.path.basename(hdf_path)}")
    return None, None


# =============================================================================
# V003 GRANULE PROCESSOR
# =============================================================================

def _extract_v003_granule_id(hdf_path):
    """Extract a granule ID from a V003 HDF filename."""
    fname = os.path.splitext(os.path.basename(hdf_path))[0]
    # Pattern: AST_XX_NNNMMDDYYYYHHMMSS_YYYYMMDDHHMMSS_NNNNN
    # Return the full base name as the ID
    return fname


def _process_v003_granule(hdf_path, product_key, out_tif_dir, out_envi_dir,
                          write_geotiff, write_envi, overwrite=True):
    """
    Process a single V003 HDF file: read bands, apply scaling, geolocate,
    and write outputs.
    """
    # Read science bands
    arrays, descriptions, wavelengths, shape = _read_v003_hdf(hdf_path, product_key)

    if not arrays:
        print(f"  [skip] No science bands read from {os.path.basename(hdf_path)}")
        return 0

    # Get geolocation
    transform, crs = _get_v003_geolocation(hdf_path, product_key, shape)

    rows, cols = shape
    files_written = 0
    gid = _extract_v003_granule_id(hdf_path)

    # Build reference metadata for GeoTIFF writing
    ref_meta = {
        'driver': 'GTiff',
        'height': rows,
        'width':  cols,
        'count':  len(arrays),
        'dtype':  'float32',
    }
    if crs:
        ref_meta['crs'] = crs
    if transform:
        ref_meta['transform'] = transform

    crs_wkt = crs.to_wkt() if crs else None

    # Write GeoTIFF
    if write_geotiff and out_tif_dir:
        out_path = os.path.join(out_tif_dir, f"{gid}.tif")
        result = _write_geotiff(arrays, out_path, ref_meta, descriptions,
                                overwrite=overwrite)
        print(f"  [GeoTIFF] {os.path.basename(result)}")
        files_written += 1

    # Write ENVI
    if write_envi and out_envi_dir:
        out_base = os.path.join(out_envi_dir, gid)
        result = _write_envi(arrays, out_base, descriptions,
                             wavelengths, crs_wkt, transform,
                             overwrite=overwrite)
        print(f"  [ENVI]    {os.path.basename(result)}")
        files_written += 1

    return files_written


# =============================================================================
# V003 MERGED VNIR+SWIR PRODUCT (AST_07M)
# =============================================================================

def _process_v003_merged(vnir_path, swir_path, out_tif_dir, out_envi_dir,
                          write_geotiff, write_envi, overwrite=True):
    """
    Create a merged 9-band AST_07M product from a VNIR + SWIR pair.
    VNIR bands are downsampled to SWIR resolution via 2×2 block averaging.
    """
    print(f"  Merging VNIR + SWIR into AST_07M...")

    # Read SWIR bands (these define the target resolution)
    swir_arrays, swir_desc, swir_wl, swir_shape = _read_v003_hdf(
        swir_path, 'AST_07')
    if not swir_arrays:
        print(f"  [skip] No SWIR bands read")
        return 0

    # Read VNIR bands
    vnir_arrays, vnir_desc, vnir_wl, vnir_shape = _read_v003_hdf(
        vnir_path, 'AST_07')
    if not vnir_arrays:
        print(f"  [skip] No VNIR bands read")
        return 0

    # Downsample VNIR to SWIR resolution (2×2 block average)
    vnir_downsampled = []
    for arr in vnir_arrays:
        ds = _downsample_2x2(arr)
        # Ensure exact shape match with SWIR (trim if off by 1 pixel)
        target_rows, target_cols = swir_shape
        ds = ds[:target_rows, :target_cols]
        vnir_downsampled.append(ds)

    # Stack: VNIR bands first (1, 2, 3N), then SWIR bands (4–9)
    merged_arrays = vnir_downsampled + swir_arrays
    merged_desc = vnir_desc + swir_desc
    merged_wl = vnir_wl + swir_wl

    # Geolocation from the SWIR file (matches target resolution)
    transform, crs = _get_v003_geolocation(swir_path, 'AST_07', swir_shape)

    rows, cols = swir_shape
    files_written = 0

    # Build granule ID from the VNIR filename, replacing AST_07 with AST_07M
    vnir_gid = _extract_v003_granule_id(vnir_path)
    merged_gid = re.sub(r'^AST_07(?:XT)?', 'AST_07M', vnir_gid,
                        flags=re.IGNORECASE)

    ref_meta = {
        'driver': 'GTiff',
        'height': rows,
        'width':  cols,
        'count':  len(merged_arrays),
        'dtype':  'float32',
    }
    if crs:
        ref_meta['crs'] = crs
    if transform:
        ref_meta['transform'] = transform

    crs_wkt = crs.to_wkt() if crs else None

    # Write GeoTIFF
    if write_geotiff and out_tif_dir:
        out_path = os.path.join(out_tif_dir, f"{merged_gid}.tif")
        result = _write_geotiff(merged_arrays, out_path, ref_meta, merged_desc,
                                overwrite=overwrite)
        print(f"  [GeoTIFF] {os.path.basename(result)}  (9 bands, merged)")
        files_written += 1

    # Write ENVI
    if write_envi and out_envi_dir:
        out_base = os.path.join(out_envi_dir, merged_gid)
        result = _write_envi(merged_arrays, out_base, merged_desc,
                             merged_wl, crs_wkt, transform,
                             overwrite=overwrite)
        print(f"  [ENVI]    {os.path.basename(result)}  (9 bands, merged)")
        files_written += 1

    return files_written


# =============================================================================
# V003 PUBLIC ENTRY POINT
# =============================================================================

def process_directory_v003(input_dir, output_dir=None,
                           aster_product='AST_05',
                           write_geotiff=True, write_envi=True,
                           overwrite=True):
    """
    Convert all ASTER V003 HDF granules in input_dir to analysis-ready files.

    Parameters
    ----------
    input_dir    : str   Path to folder containing .hdf files (flat or subfolders).
    output_dir   : str   Output root.  Defaults to <input_dir>/converted/ASTER_V003.
    aster_product: str   One of: 'AST_05', 'AST_07', 'AST_07M', 'AST_08', 'AST_09T'.
                         'AST_07M' processes all AST_07/07XT individually AND creates
                         a merged 9-band product for matched VNIR+SWIR pairs.
                         If 'auto', auto-detects from filenames.
    write_geotiff: bool  Write GeoTIFF output.
    write_envi   : bool  Write ENVI output.
    overwrite    : bool  If True, overwrite existing files. If False, auto-increment.
    """
    if not _HAS_PYHDF:
        print("ERROR: pyhdf is required for V003 HDF files.  "
              "Install with: conda install -c conda-forge pyhdf")
        return

    # Normalise product key
    merged_mode = False
    if aster_product in ('AST_07M',):
        product_key = 'AST_07'
        display_name = 'AST_07M (Merged VNIR+SWIR)'
        merged_mode = True
    elif aster_product in ('AST_07', 'AST_07XT', 'AST_07X'):
        product_key = 'AST_07'
        display_name = 'AST_07'
    else:
        product_key = aster_product
        display_name = aster_product

    auto_detect = (product_key == 'auto')

    if not auto_detect and product_key not in V003_PRODUCTS:
        print(f"ERROR: Unknown V003 product '{aster_product}'.  "
              f"Choose from: AST_05, AST_07, AST_07M, AST_08, AST_09T")
        return

    if not write_geotiff and not write_envi:
        print("ERROR: At least one of write_geotiff or write_envi must be True.")
        return

    # Set up output directories
    if output_dir is None:
        output_dir = os.path.join(input_dir, 'converted', 'ASTER_V003')

    out_tif_dir  = os.path.join(output_dir, 'GeoTIFF') if write_geotiff else None
    out_envi_dir = os.path.join(output_dir, 'ENVI')    if write_envi    else None

    for d in filter(None, (out_tif_dir, out_envi_dir)):
        os.makedirs(d, exist_ok=True)

    print(f"ASTER V003 HDF Converter  |  product: "
          f"{'auto-detect' if auto_detect else display_name}")
    print(f"Input  : {input_dir}")
    print(f"Output : {output_dir}")
    print(f"Format : {'GeoTIFF ' if write_geotiff else ''}{'ENVI' if write_envi else ''}")
    print("-" * 60)

    # ------------------------------------------------------------------
    # Discover HDF files — scan subfolders first, then flat folder
    # ------------------------------------------------------------------
    hdf_files = []

    # Check subfolders
    for entry in sorted(os.listdir(input_dir)):
        sf = os.path.join(input_dir, entry)
        if os.path.isdir(sf):
            for f in glob.glob(os.path.join(sf, '*.hdf')):
                hdf_files.append(f)
            for f in glob.glob(os.path.join(sf, '*.HDF')):
                hdf_files.append(f)

    # Flat folder
    for f in glob.glob(os.path.join(input_dir, '*.hdf')):
        if f not in hdf_files:
            hdf_files.append(f)
    for f in glob.glob(os.path.join(input_dir, '*.HDF')):
        if f not in hdf_files:
            hdf_files.append(f)

    if not hdf_files:
        print("ERROR: No .hdf files found in input directory.")
        return

    # Filter by product type
    granules = []
    for hdf_path in sorted(hdf_files):
        detected = _detect_v003_product(hdf_path)
        if auto_detect:
            if detected and detected in V003_PRODUCTS:
                granules.append((hdf_path, detected))
        else:
            if detected == product_key:
                granules.append((hdf_path, product_key))

    if not granules:
        print(f"No {'V003 HDF' if auto_detect else display_name} files found.")
        return

    print(f"Found {len(granules)} V003 HDF granule(s)")
    print("-" * 60)

    total_files = 0
    for i, (hdf_path, pkey) in enumerate(granules, 1):
        fname = os.path.basename(hdf_path)
        pname = V003_PRODUCTS[pkey]['name']
        print(f"\n[{i}/{len(granules)}] {fname}")
        print(f"  Product: {pname}")
        n = _process_v003_granule(hdf_path, pkey,
                                  out_tif_dir, out_envi_dir,
                                  write_geotiff, write_envi,
                                  overwrite=overwrite)
        total_files += n

    # --- AST_07M: match VNIR+SWIR pairs and create merged products ---
    if merged_mode and product_key == 'AST_07':
        # Separate VNIR and SWIR files by reading their SDS contents
        vnir_files = {}  # acquisition_id -> path
        swir_files = {}

        for hdf_path, _ in granules:
            acq_id = _extract_v003_acquisition_id(hdf_path)
            if acq_id is None:
                continue
            ftype = _detect_v003_07_type(hdf_path)
            if ftype == 'vnir':
                vnir_files[acq_id] = hdf_path
            elif ftype == 'swir':
                swir_files[acq_id] = hdf_path
            elif ftype == 'both':
                # File has both VNIR and SWIR — use as both
                vnir_files[acq_id] = hdf_path
                swir_files[acq_id] = hdf_path

        # Find matched pairs
        matched_ids = set(vnir_files.keys()) & set(swir_files.keys())

        if matched_ids:
            print(f"\n{'=' * 60}")
            print(f"Creating {len(matched_ids)} merged AST_07M product(s)...")
            print("=" * 60)

            for j, acq_id in enumerate(sorted(matched_ids), 1):
                vnir_path = vnir_files[acq_id]
                swir_path = swir_files[acq_id]
                print(f"\n[Merge {j}/{len(matched_ids)}]")
                print(f"  VNIR: {os.path.basename(vnir_path)}")
                print(f"  SWIR: {os.path.basename(swir_path)}")
                n = _process_v003_merged(vnir_path, swir_path,
                                          out_tif_dir, out_envi_dir,
                                          write_geotiff, write_envi,
                                          overwrite=overwrite)
                total_files += n
        else:
            unmatched_vnir = set(vnir_files.keys()) - set(swir_files.keys())
            unmatched_swir = set(swir_files.keys()) - set(vnir_files.keys())
            if unmatched_vnir or unmatched_swir:
                print(f"\n  No matched VNIR+SWIR pairs found for merging.")
                if unmatched_vnir:
                    print(f"  VNIR-only files: {len(unmatched_vnir)}")
                if unmatched_swir:
                    print(f"  SWIR-only files: {len(unmatched_swir)}")

    print("-" * 60)
    print(f"Done.  {total_files} file(s) written to: {output_dir}")




if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("ASTER Converter \u2014 part of RSDTK")
        print("=" * 50)
        print(f"\nUsage: python {os.path.basename(__file__)} <input_dir> [output_dir] [options]")
        print("\nProduct:")
        print("  --product <n>   AST_05 | AST_07 | AST_08 | AST_09T | AST_L1T  (default: AST_05)")
        print("\nData version:")
        print("  --v003           Process V003 HDF files (requires pyhdf)")
        print("  --v004           Process V004 GeoTIFF files (default)")
        print("\nOutput format:")
        print("  --geotiff          GeoTIFF only")
        print("  --envi             ENVI only")
        print("  (default: both)")
        print("\nOverwrite:")
        print("  --overwrite        Overwrite existing output files (default)")
        print("  --no-overwrite     Auto-increment filenames to avoid overwriting")
        print("\nExamples:")
        print(f"  python {os.path.basename(__file__)} D:\\Data\\ASTER")
        print(f"  python {os.path.basename(__file__)} D:\\Data\\ASTER --product AST_09T --v003")
        print(f"  python {os.path.basename(__file__)} D:\\Data\\ASTER --product AST_09T --envi --no-overwrite")
        sys.exit(0)

    args      = sys.argv[1:]
    positional = []
    product    = 'AST_05'
    write_tif  = True
    write_envi = True
    use_v003   = False
    overwrite  = True

    i = 0
    while i < len(args):
        a = args[i].lower()
        if a == '--product' and i + 1 < len(args):
            product = args[i + 1].upper()
            i += 2
        elif a == '--geotiff':
            write_tif  = True
            write_envi = False
            i += 1
        elif a == '--envi':
            write_tif  = False
            write_envi = True
            i += 1
        elif a == '--v003':
            use_v003 = True
            i += 1
        elif a == '--v004':
            use_v003 = False
            i += 1
        elif a == '--overwrite':
            overwrite = True
            i += 1
        elif a == '--no-overwrite':
            overwrite = False
            i += 1
        else:
            positional.append(args[i])
            i += 1

    input_dir  = positional[0]
    output_dir = positional[1] if len(positional) > 1 else None

    if not os.path.isdir(input_dir):
        print(f"ERROR: Input directory does not exist: {input_dir}")
        sys.exit(1)

    if use_v003:
        process_directory_v003(input_dir, output_dir,
                               aster_product=product,
                               write_geotiff=write_tif,
                               write_envi=write_envi,
                               overwrite=overwrite)
    else:
        process_directory(input_dir, output_dir,
                          aster_product=product,
                          write_geotiff=write_tif,
                          write_envi=write_envi,
                          overwrite=overwrite)