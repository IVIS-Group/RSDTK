"""
Spectral Resampling and Band Selection Module
Universal tool for hyperspectral data (EMIT, MASTER, HyTES, AVIRIS-NG, etc.)

Features:
  - Wavelength region selection (VIS, VNIR, SWIR, full)
  - Match to multispectral sensors (Landsat, Sentinel-2, ASTER, MODIS)
  - Every-Nth band decimation
  - Custom wavelength range selection
  - Water vapor absorption band exclusion
  - Gaussian spectral convolution to target sensor bands
  - RSR file loading for accurate convolution (when available)

All wavelengths in nanometers unless otherwise noted.
"""

import os
import sys
import numpy as np
from collections import OrderedDict


# ============================================================================
# Water vapor absorption regions to auto-exclude (nm)
# ============================================================================

WATER_VAPOR_REGIONS = [
    (1334, 1431),   # Strong H2O absorption
    (1788, 1980),   # Strong H2O absorption
]


# ============================================================================
# Wavelength region presets (nm)
# ============================================================================

WAVELENGTH_REGIONS = {
    'vis':  (380, 700),     # Visible only
    'vnir': (380, 1000),    # Visible + Near-IR
    'nir':  (700, 1000),    # Near-IR only
    'swir': (1000, 2500),   # Shortwave IR only
    'full': (380, 2500),    # Full EMIT range
}


# ============================================================================
# Multispectral sensor band definitions
# Center wavelength (nm), bandwidth/FWHM (nm), band name
# Only VSWIR bands included (no TIR since EMIT doesn't cover TIR)
# ============================================================================

SENSOR_BANDS = {
    'landsat89': OrderedDict([
        ('B1_Coastal',  {'center': 443, 'fwhm': 20,  'range': (433, 453)}),
        ('B2_Blue',     {'center': 482, 'fwhm': 65,  'range': (450, 515)}),
        ('B3_Green',    {'center': 562, 'fwhm': 75,  'range': (525, 600)}),
        ('B4_Red',      {'center': 655, 'fwhm': 50,  'range': (630, 680)}),
        ('B5_NIR',      {'center': 865, 'fwhm': 30,  'range': (845, 885)}),
        ('B6_SWIR1',    {'center': 1609, 'fwhm': 80, 'range': (1566, 1651)}),
        ('B7_SWIR2',    {'center': 2201, 'fwhm': 180,'range': (2107, 2294)}),
    ]),

    'sentinel2': OrderedDict([
        ('B1_Coastal',  {'center': 443,  'fwhm': 20,  'range': (433, 453)}),
        ('B2_Blue',     {'center': 490,  'fwhm': 65,  'range': (458, 523)}),
        ('B3_Green',    {'center': 560,  'fwhm': 35,  'range': (543, 578)}),
        ('B4_Red',      {'center': 665,  'fwhm': 30,  'range': (650, 680)}),
        ('B5_RedEdge1', {'center': 705,  'fwhm': 15,  'range': (698, 713)}),
        ('B6_RedEdge2', {'center': 740,  'fwhm': 15,  'range': (733, 748)}),
        ('B7_RedEdge3', {'center': 783,  'fwhm': 20,  'range': (773, 793)}),
        ('B8_NIR',      {'center': 842,  'fwhm': 115, 'range': (785, 900)}),
        ('B8A_NIRn',    {'center': 865,  'fwhm': 20,  'range': (855, 875)}),
        ('B9_WaterVap', {'center': 945,  'fwhm': 20,  'range': (935, 955)}),
        ('B11_SWIR1',   {'center': 1610, 'fwhm': 90,  'range': (1565, 1655)}),
        ('B12_SWIR2',   {'center': 2190, 'fwhm': 180, 'range': (2100, 2280)}),
    ]),

    'aster': OrderedDict([
        ('B1_Green',  {'center': 556,  'fwhm': 80,  'range': (520, 600)}),
        ('B2_Red',    {'center': 661,  'fwhm': 60,  'range': (630, 690)}),
        ('B3N_NIR',   {'center': 807,  'fwhm': 80,  'range': (760, 860)}),
        ('B4_SWIR1',  {'center': 1656, 'fwhm': 100, 'range': (1600, 1700)}),
        ('B5_SWIR2',  {'center': 2167, 'fwhm': 40,  'range': (2145, 2185)}),
        ('B6_SWIR3',  {'center': 2209, 'fwhm': 40,  'range': (2185, 2225)}),
        ('B7_SWIR4',  {'center': 2262, 'fwhm': 50,  'range': (2235, 2285)}),
        ('B8_SWIR5',  {'center': 2336, 'fwhm': 70,  'range': (2295, 2365)}),
        ('B9_SWIR6',  {'center': 2400, 'fwhm': 70,  'range': (2360, 2430)}),
    ]),

    'modis_land': OrderedDict([
        ('B1_Red',    {'center': 645,  'fwhm': 50,  'range': (620, 670)}),
        ('B2_NIR',    {'center': 858,  'fwhm': 35,  'range': (841, 876)}),
        ('B3_Blue',   {'center': 469,  'fwhm': 20,  'range': (459, 479)}),
        ('B4_Green',  {'center': 555,  'fwhm': 20,  'range': (545, 565)}),
        ('B5_SWIR1',  {'center': 1240, 'fwhm': 20,  'range': (1230, 1250)}),
        ('B6_SWIR2',  {'center': 1640, 'fwhm': 24,  'range': (1628, 1652)}),
        ('B7_SWIR3',  {'center': 2130, 'fwhm': 50,  'range': (2105, 2155)}),
    ]),
}

# Sensor name aliases for user convenience
SENSOR_ALIASES = {
    'landsat': 'landsat89',
    'landsat8': 'landsat89',
    'landsat9': 'landsat89',
    'l8': 'landsat89',
    'l9': 'landsat89',
    'oli': 'landsat89',
    's2': 'sentinel2',
    'sentinel': 'sentinel2',
    'msi': 'sentinel2',
    'modis': 'modis_land',
}


# ============================================================================
# Band selection functions
# ============================================================================

def get_good_band_mask(good_wavelengths, wavelengths, exclude_water=True):
    """
    Build a boolean mask for bands to include.

    Parameters:
        good_wavelengths: array of 0/1 flags from EMIT
        wavelengths: array of center wavelengths (nm)
        exclude_water: also exclude water vapor absorption regions

    Returns:
        boolean array (True = include)
    """
    mask = good_wavelengths.astype(bool)

    if exclude_water:
        for wv_min, wv_max in WATER_VAPOR_REGIONS:
            water_region = (wavelengths >= wv_min) & (wavelengths <= wv_max)
            mask = mask & ~water_region
            n_excluded = water_region.sum()
            if n_excluded > 0:
                print(f"  Excluding {n_excluded} bands in water vapor region "
                      f"{wv_min}-{wv_max} nm")

    return mask


def select_wavelength_region(wavelengths, region_name):
    """
    Select bands within a named wavelength region.

    Parameters:
        wavelengths: array of center wavelengths (nm)
        region_name: key from WAVELENGTH_REGIONS dict

    Returns:
        boolean mask array
    """
    region_name = region_name.lower()
    if region_name not in WAVELENGTH_REGIONS:
        raise ValueError(f"Unknown region '{region_name}'. "
                         f"Available: {list(WAVELENGTH_REGIONS.keys())}")

    wl_min, wl_max = WAVELENGTH_REGIONS[region_name]
    mask = (wavelengths >= wl_min) & (wavelengths <= wl_max)
    print(f"  Region '{region_name}': {wl_min}-{wl_max} nm, "
          f"{mask.sum()} bands selected")
    return mask


def select_wavelength_range(wavelengths, ranges_str):
    """
    Select bands within user-specified wavelength ranges.

    Parameters:
        wavelengths: array of center wavelengths (nm)
        ranges_str: comma-separated ranges like "450-900,2000-2400"

    Returns:
        boolean mask array
    """
    mask = np.zeros(len(wavelengths), dtype=bool)

    for range_part in ranges_str.split(','):
        range_part = range_part.strip()
        if '-' in range_part:
            parts = range_part.split('-')
            wl_min = float(parts[0])
            wl_max = float(parts[1])
            region = (wavelengths >= wl_min) & (wavelengths <= wl_max)
            mask = mask | region
            print(f"  Custom range {wl_min}-{wl_max} nm: "
                  f"{region.sum()} bands")

    print(f"  Total bands selected: {mask.sum()}")
    return mask


def select_every_nth(wavelengths, n):
    """
    Select every Nth band.

    Parameters:
        wavelengths: array of center wavelengths (nm)
        n: step size

    Returns:
        boolean mask array
    """
    mask = np.zeros(len(wavelengths), dtype=bool)
    mask[::n] = True
    print(f"  Every {n}th band: {mask.sum()} of {len(wavelengths)} bands")
    return mask


# ============================================================================
# Spectral resampling to target sensor
# ============================================================================

def gaussian_rsr(wavelengths, center, fwhm):
    """
    Generate a Gaussian relative spectral response function.

    Parameters:
        wavelengths: array of source wavelengths (nm)
        center: center wavelength of target band (nm)
        fwhm: full width at half maximum (nm)

    Returns:
        array of weights (normalized to sum to 1)
    """
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    weights = np.exp(-0.5 * ((wavelengths - center) / sigma) ** 2)

    # Zero out weights far from center (beyond 2x FWHM)
    weights[(wavelengths < center - fwhm) | (wavelengths > center + fwhm)] = 0

    weight_sum = weights.sum()
    if weight_sum > 0:
        weights /= weight_sum
    return weights


def load_rsr_file(filepath):
    """
    Load a relative spectral response file.
    Expected format: CSV or whitespace-delimited with columns:
      wavelength_nm, response

    Parameters:
        filepath: path to RSR file

    Returns:
        (wavelengths, response) tuple of arrays
    """
    data = np.loadtxt(filepath, delimiter=None, comments='#')
    if data.shape[1] >= 2:
        return data[:, 0], data[:, 1]
    else:
        raise ValueError(f"RSR file must have at least 2 columns: {filepath}")


def resample_to_sensor(data_cube, src_wavelengths, target_sensor,
                       rsr_dir=None):
    """
    Resample a hyperspectral cube to match a target multispectral sensor.

    Parameters:
        data_cube: array with shape (lines, samples, bands) or (bands, lines, samples)
        src_wavelengths: array of source wavelengths (nm)
        target_sensor: sensor name string (e.g., 'landsat89', 'sentinel2')
        rsr_dir: optional directory containing RSR files named
                 {sensor}_{band_name}.csv

    Returns:
        resampled_cube: 3D array with target sensor bands
        band_names: list of target band names
        band_wavelengths: list of center wavelengths (nm)
    """
    # Resolve sensor alias
    sensor_key = SENSOR_ALIASES.get(target_sensor.lower(), target_sensor.lower())
    if sensor_key not in SENSOR_BANDS:
        raise ValueError(f"Unknown sensor '{target_sensor}'. "
                         f"Available: {list(SENSOR_BANDS.keys())} "
                         f"(aliases: {list(SENSOR_ALIASES.keys())})")

    bands = SENSOR_BANDS[sensor_key]

    # Determine data layout
    band_last = (data_cube.shape[2] == len(src_wavelengths))
    if not band_last:
        # Assume (bands, lines, samples) - transpose to (lines, samples, bands)
        data_cube = np.transpose(data_cube, (1, 2, 0))

    n_lines, n_samples, n_src_bands = data_cube.shape
    n_target = len(bands)

    resampled = np.full((n_lines, n_samples, n_target), np.nan, dtype=np.float32)
    band_names = []
    band_wavelengths = []

    print(f"  Resampling to {sensor_key} ({n_target} bands):")

    for i, (bname, binfo) in enumerate(bands.items()):
        center = binfo['center']
        fwhm = binfo['fwhm']

        # Try to load RSR file if directory provided
        rsr_weights = None
        if rsr_dir:
            rsr_file = os.path.join(rsr_dir, f"{sensor_key}_{bname}.csv")
            if os.path.exists(rsr_file):
                rsr_wl, rsr_resp = load_rsr_file(rsr_file)
                # Interpolate RSR to source wavelengths
                rsr_weights = np.interp(src_wavelengths, rsr_wl, rsr_resp,
                                        left=0, right=0)
                rsr_sum = rsr_weights.sum()
                if rsr_sum > 0:
                    rsr_weights /= rsr_sum
                print(f"    {bname}: {center} nm (RSR from file)")

        # Fall back to Gaussian if no RSR file
        if rsr_weights is None:
            rsr_weights = gaussian_rsr(src_wavelengths, center, fwhm)
            contributing = (rsr_weights > 0).sum()
            print(f"    {bname}: {center} nm ± {fwhm/2:.0f} nm "
                  f"({contributing} source bands)")

        # Check if any source bands contribute
        if rsr_weights.sum() == 0:
            print(f"    WARNING: No source bands cover {bname} ({center} nm)")
            band_names.append(bname)
            band_wavelengths.append(center)
            continue

        # Apply weighted sum
        resampled[:, :, i] = np.nansum(
            data_cube * rsr_weights[np.newaxis, np.newaxis, :], axis=2
        )

        band_names.append(bname)
        band_wavelengths.append(center)

    return resampled, band_names, np.array(band_wavelengths)


def nearest_band_match(src_wavelengths, target_sensor):
    """
    Find the nearest source band index for each target sensor band.
    Used as a simple alternative to Gaussian convolution.

    Parameters:
        src_wavelengths: array of source wavelengths (nm)
        target_sensor: sensor name string

    Returns:
        indices: array of source band indices
        band_names: list of target band names
        band_wavelengths: array of target center wavelengths
    """
    sensor_key = SENSOR_ALIASES.get(target_sensor.lower(), target_sensor.lower())
    if sensor_key not in SENSOR_BANDS:
        raise ValueError(f"Unknown sensor '{target_sensor}'")

    bands = SENSOR_BANDS[sensor_key]
    indices = []
    band_names = []
    band_wavelengths = []

    print(f"  Nearest-band matching to {sensor_key}:")
    for bname, binfo in bands.items():
        center = binfo['center']
        idx = np.argmin(np.abs(src_wavelengths - center))
        actual_wl = src_wavelengths[idx]
        offset = actual_wl - center
        print(f"    {bname}: target {center} nm -> source {actual_wl:.1f} nm "
              f"(offset {offset:+.1f} nm)")
        indices.append(idx)
        band_names.append(bname)
        band_wavelengths.append(center)

    return np.array(indices), band_names, np.array(band_wavelengths)


# ============================================================================
# Command-line argument parsing
# ============================================================================

def parse_spectral_args(args):
    """
    Parse spectral resampling command line arguments.

    Recognized arguments:
        --region <name>          Select wavelength region (vis/vnir/nir/swir/full)
        --bands <ranges>         Custom wavelength ranges (e.g., 450-900,2000-2400)
        --match-sensor <name>    Resample to sensor (landsat/sentinel2/aster/modis)
        --nearest                Use nearest band instead of Gaussian convolution
        --every-nth <n>          Select every Nth band
        --keep-water             Don't exclude water vapor bands
        --rsr-dir <path>         Directory with RSR files for convolution

    Returns:
        dict with parsed options
    """
    result = {
        'mode': 'default',        # default/region/bands/sensor/nth
        'region': None,
        'ranges': None,
        'sensor': None,
        'nearest': False,
        'every_nth': None,
        'exclude_water': True,
        'rsr_dir': None,
    }

    i = 0
    while i < len(args):
        arg = args[i].lower()

        if arg == '--region':
            if i + 1 >= len(args):
                raise ValueError("--region requires a name (vis/vnir/nir/swir/full)")
            result['mode'] = 'region'
            result['region'] = args[i + 1].lower()
            i += 2

        elif arg == '--bands':
            if i + 1 >= len(args):
                raise ValueError("--bands requires wavelength ranges")
            result['mode'] = 'bands'
            result['ranges'] = args[i + 1]
            i += 2

        elif arg == '--match-sensor':
            if i + 1 >= len(args):
                raise ValueError("--match-sensor requires a sensor name")
            result['mode'] = 'sensor'
            result['sensor'] = args[i + 1]
            i += 2

        elif arg == '--nearest':
            result['nearest'] = True
            i += 1

        elif arg == '--every-nth':
            if i + 1 >= len(args):
                raise ValueError("--every-nth requires a number")
            result['mode'] = 'nth'
            result['every_nth'] = int(args[i + 1])
            i += 2

        elif arg == '--keep-water':
            result['exclude_water'] = False
            i += 1

        elif arg == '--rsr-dir':
            if i + 1 >= len(args):
                raise ValueError("--rsr-dir requires a directory path")
            result['rsr_dir'] = args[i + 1]
            i += 2

        else:
            i += 1

    return result


def apply_spectral_selection(wavelengths, good_wavelengths, spectral_opts):
    """
    Apply spectral selection options to get a band mask.

    Parameters:
        wavelengths: array of center wavelengths (nm)
        good_wavelengths: array of 0/1 quality flags
        spectral_opts: dict from parse_spectral_args()

    Returns:
        band_mask: boolean array of bands to include
        resample_info: dict with sensor resampling info (or None)
    """
    exclude_water = spectral_opts.get('exclude_water', True)
    resample_info = None

    # Start with good wavelength mask
    base_mask = get_good_band_mask(good_wavelengths, wavelengths, exclude_water)

    mode = spectral_opts.get('mode', 'default')

    if mode == 'default':
        return base_mask, None

    elif mode == 'region':
        region_mask = select_wavelength_region(wavelengths, spectral_opts['region'])
        return base_mask & region_mask, None

    elif mode == 'bands':
        range_mask = select_wavelength_range(wavelengths, spectral_opts['ranges'])
        return base_mask & range_mask, None

    elif mode == 'nth':
        nth_mask = select_every_nth(wavelengths, spectral_opts['every_nth'])
        return base_mask & nth_mask, None

    elif mode == 'sensor':
        # For sensor matching, return all good bands (resampling happens later)
        resample_info = {
            'sensor': spectral_opts['sensor'],
            'nearest': spectral_opts.get('nearest', False),
            'rsr_dir': spectral_opts.get('rsr_dir'),
        }
        return base_mask, resample_info

    return base_mask, None


# ============================================================================
# Utility: list available sensors
# ============================================================================

def list_sensors():
    """Print available sensor definitions."""
    print("Available sensor definitions for --match-sensor:")
    print("-" * 60)
    for sensor_key, bands in SENSOR_BANDS.items():
        aliases = [k for k, v in SENSOR_ALIASES.items() if v == sensor_key]
        alias_str = f" (aliases: {', '.join(aliases)})" if aliases else ""
        print(f"\n  {sensor_key}{alias_str}:")
        for bname, binfo in bands.items():
            print(f"    {bname}: {binfo['center']} nm "
                  f"(FWHM={binfo['fwhm']} nm, "
                  f"range={binfo['range'][0]}-{binfo['range'][1]} nm)")


# ============================================================================
# Main - standalone testing/info
# ============================================================================

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] == '--help':
        print("Spectral Resampling and Band Selection Module")
        print("=" * 50)
        print(f"\nUsage: python {os.path.basename(__file__)} --list-sensors")
        print(f"\nSpectral options (used with EMIT/MASTER/HyTES converters):")
        print(f"  --region <name>        Wavelength region: vis, vnir, nir, swir, full")
        print(f"  --bands <ranges>       Custom ranges: e.g., 450-900,2000-2400")
        print(f"  --match-sensor <name>  Resample to sensor bands")
        print(f"  --nearest              Use nearest band (not Gaussian convolution)")
        print(f"  --every-nth <n>        Keep every Nth band")
        print(f"  --keep-water           Don't auto-exclude water vapor bands")
        print(f"  --rsr-dir <path>       Directory with RSR CSV files")
        print(f"\nWater vapor regions auto-excluded:")
        for wv_min, wv_max in WATER_VAPOR_REGIONS:
            print(f"  {wv_min}-{wv_max} nm")
        print(f"\nWavelength regions:")
        for name, (wl_min, wl_max) in WAVELENGTH_REGIONS.items():
            print(f"  {name}: {wl_min}-{wl_max} nm")
        sys.exit(0)

    if sys.argv[1] == '--list-sensors':
        list_sensors()