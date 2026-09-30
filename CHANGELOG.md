# Changelog

Notable changes to RSDTK. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and from v1.3 onward
versioning follows [Semantic Versioning](https://semver.org/):

- **Patch** (1.2.x) — bug fixes only, no behavior change
- **Minor** (1.x.0) — new sensors or features, existing outputs unchanged
- **Major** (x.0.0) — changes that alter existing output or break a workflow

> **On development history.** RSDTK grew out of an ASTER-only toolkit begun in
> February 2026, became a unified multi-sensor application in March 2026, and
> was released publicly in October 2026. Versions 1.0 to 1.2 therefore predate
> the public repository, and the entries below are reconstructed from
> development records and the user guide's version notes rather than from commit
> history. Dates are approximate.

---

## [Unreleased]

Planned work is tracked in [ROADMAP.md](ROADMAP.md).

### Fixed

- **MODIS ENVI output could not be georeferenced correctly in ENVI.** The header
  wrote a `coordinate system string` containing a GDAL-normalized WKT in which
  the EPSG authority code contradicted the inline sphere radius: the spheroid was
  declared as `SPHEROID["WGS 84", 6371007.181, 0, AUTHORITY["EPSG","7030"]]`,
  where 6371007.181 m is the correct MODIS authalic sphere but EPSG:7030 is the
  actual WGS 84 ellipsoid (6378137 m, 1/298.257). ENVI resolves the authority
  code in preference to the inline value and so inverted the sinusoidal
  projection on an ellipsoid rather than a sphere, reporting geographic
  coordinates that were out by thousands of degrees.

  The header now writes an explicit `projection info` line giving the projection
  type, sphere radius and central meridian, adds `units=Meters` to `map info`,
  and omits the conflicting WKT entirely. Together `map info` and
  `projection info` fully define the grid for ENVI.

  GeoTIFF output was never affected, since ArcGIS and QGIS tolerate the
  contradictory authority code and use the inline radius.

### Changed

- **Source tree reconstructed after storage corruption.** Twelve source files
  were corrupted by a drive failure and subsequent filesystem recovery. Each was
  restored from the bytecode bundled in the v1.2 installer and verified against
  it: definition inventories match exactly across all twelve, as do every
  `process_directory()` signature. Logic and constants are faithful to the
  shipped v1.2 build; comments and docstrings in the affected files were lost to
  compilation and have been rewritten where present. Affected: `ASTER_Converter`,
  `AVIRIS_Converter`, `Landsat_L1_Converter`, `Landsat_L2_Converter`,
  `MODIS_Converter`, `Sentinel2_L2A_Converter`, `Sentinel3_Converter`,
  `VIIRS_Converter`, `spectral_tools`, `subset_tool`, `runtime_hook`, and
  `build_rsdtk.bat`.

---

## [1.2.0] — May 2026

The release that took RSDTK from a handful of sensors to 14 sensor families and
30+ product types.

### Added

- **MODIS support.** MOD09GA/MYD09GA surface reflectance (bands 1–7, 500 m) and
  MOD11A1/MYD11A1 LST and emissivity (1 km). Native sinusoidal projection
  retained, using a WKT-based CRS to sidestep PROJ database issues.
- **Sentinel-3 support.** SLSTR L1B RBT, SLSTR L2 LST, OLCI L1B EFR (21 bands),
  and Synergy L2 SYN surface reflectance. KD-tree reprojection from swath to
  geographic WGS-84, with a `product_filter` parameter so selecting one product
  processes only matching `.SEN3` folders.
- **AVIRIS-NG and Classic AVIRIS support.** L1B radiance and L2 reflectance from
  orthorectified ENVI binaries.
- **AVIRIS-3/5 L2A support.** Orthorectified surface reflectance from
  EarthData/ORNL DAAC NetCDF, already on a UTM grid so no GLT step needed. Read
  from NetCDF groups, with pyproj-based bbox cropping converting lon/lat to UTM.
- **AVIRIS-5 detection**, alongside AVIRIS-3, from the filename prefix.
- **Landsat, all missions.** Extended from Landsat 8/9 to the full archive:
  OLI/TIRS (L8–9), ETM+ (L7), TM (L4–5), and MSS (L1–5), auto-detected from MTL
  metadata.

### Changed

- **GUI restructured** from 15 sensor cards to 10, with product selection moved
  into dropdowns. MODIS and Sentinel-3 cards added; the Landsat card renamed and
  given file-prefix descriptions.
- **Water vapor exclusion** widened from 1334–1431 nm to 1334–1440 nm across all
  AVIRIS converters, catching edge channels that returned negative reflectance.

### Fixed

- `discover_scenes()` selecting the wrong MTL file where L1TP and L2SP products
  share a folder.
- Sentinel-3 resolution estimation producing grids that were either too fine or
  too coarse; now uses `max(pixel_spacing, sqrt(area / n_pixels))`.
- Sentinel-3 scale and offset being applied twice, since netCDF4 applies them
  automatically; resolved with `set_auto_maskandscale(False)` before reading.
- Synergy variable naming handled correctly (`SDR_Oa01` inside
  `Syn_Oa01_reflectance.nc`), and the absence of Oa13–15 and Oa19–20 accounted
  for.

---

## [1.1.0] — May 2026

### Added

- **QA file output for ASTER V004** as separate multiband `_QA` GeoTIFF and ENVI
  files, rather than appended to the science bands.
- **Reflectance and emissivity validity masking.** Pixels outside 0–1 are set to
  NaN in GeoTIFF and −9999 in ENVI.
- **"Open File Location" button**, opening Explorer at the output folder once
  processing finishes.

### Changed

- **Sensor card icons** replaced with consistent two-letter badge labels (AS,
  L1, L2, S2, E1, E2, VI, EM, AI, MA, HT, AV) on a fixed red accent badge.
  Emoji rendered at inconsistent sizes and broke text alignment.
- **ENVI invalid-pixel value** settled at −9999 rather than 0. Zero is a
  legitimate reflectance value, and in ratio calculations a zero in one band
  produces a spurious ±1 result — this surfaced as widespread NDVI values of 1.

### Removed

- **`_organize_output_files()`**, the post-processing step that moved outputs
  into `GeoTIFF/` and `ENVI/` subfolders. It was also moving original ASTER V004
  *input* TIFs into subfolders. Every converter already creates these
  subdirectories internally, so the step was both redundant and destructive.

### Fixed

- Overwrite parameter wiring across all converters. Automated patching had
  inserted `if not overwrite` checks inside functions such as
  `write_multiband()` without adding `overwrite=True` to those signatures or
  threading it through `process_directory()` → `process_scene()` →
  `write_multiband()`, so the toggle silently had no effect.
- JP2OpenJPEG driver not bundled in the standalone executable, which prevented
  Sentinel-2 processing from the installed build.
- Crash in windowed mode caused by `sys.stdout` being `None`, fixed with None
  checks in the `_Redirector` class.

---

## [1.0.0] — April 2026

First packaged release, distributed as a standalone Windows installer built with
PyInstaller and Inno Setup.

### Sensors at 1.0

- **ASTER V004** (GeoTIFF): AST_05, AST_07/07XT, AST_08, AST_09T, AST_L1T
- **ASTER V003** (HDF4): AST_05, AST_07, AST_08, AST_09T, and the merged AST_07M
  product
- **Landsat 8/9** Collection 2, Levels 1 and 2
- **Sentinel-2** L2A
- **ECOSTRESS** L1CT radiance and L2 LSTE
- **VIIRS** VNP02IMG, VNP02MOD, VNP21
- **EMIT** L2A
- **AIRS** L1B
- **MASTER** L1B and L2
- **HyTES** L1 and L2
- **AVIRIS-3** L1B

### Core capabilities

- Multiband GeoTIFF output, float32, DEFLATE compressed and tiled
- ENVI BSQ output with wavelength and FWHM metadata written into the header
- KD-tree swath reprojection with inverse distance weighting, using
  zero-initialized accumulation to avoid NaN propagation
- GLT orthorectification where the product provides a lookup table
- Per-product calibration and scaling factors applied automatically
- Fill-value handling: NaN in GeoTIFF, −9999 in ENVI
- Spatial subsetting by bounding box, shapefile, KML, GeoJSON, or auto-tiling
- Spectral subsetting: wavelength regions, custom ranges, sensor matching via
  Gaussian convolution or nearest band, and every-Nth decimation
- Overwrite toggle across all converters
- VIIRS band selection presets
- AVIRIS-3 spatial resolution resampling
- Separate output files per group for AST_09T (SRA/SIR) and AST_L1T
  (VNIR/SWIR/TIR)
- Threaded processing with live log output

---

## Pre-history

**February–March 2026.** RSDTK began as the ASTER Data Toolkit, an ASTER-only
application with three tools: multiband image creation, scaling factor
application, and format conversion. A standalone VIIRS Data Converter followed
in March 2026. Later that month the individual converters were unified behind a
single CustomTkinter interface, which became RSDTK.
