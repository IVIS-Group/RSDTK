# Third-Party Licenses

RSDTK is released under the MIT License (see [LICENSE](LICENSE)).

The Windows installer bundles the components below. Running RSDTK from source
installs these via conda or pip instead, in which case their licenses apply to
those installations in the normal way.

Generated from the build output in `dist\RSDTK\_internal` and the build environment's package
metadata, by `scripts/collect_licenses.py`.

| Component | Version | License | Project |
|---|---|---|---|
| Expat | 2.7.4 | MIT | https://libexpat.github.io/ |
| GDAL | 3.8.5 | MIT | https://gdal.org/ |
| GEOS | 3.12.2 | LGPL-2.1-only | https://libgeos.org/ |
| HDF4 | 4.2.15 | BSD-3-Clause | https://www.hdfgroup.org/ |
| LAPACK | 3.11.0 | BSD-3-Clause | https://www.netlib.org/lapack/ |
| NetCDF-C |  | BSD 3-Clause | https://www.unidata.ucar.edu/software/netcdf/ |
| OpenSSL | 3.6.1 | Apache-2.0 | https://www.openssl.org/ |
| PROJ | 9.4.1 | MIT | https://proj.org/ |
| SQLite | 3.51.1 | blessing | https://www.sqlite.org/ |
| Tcl |  | Tcl/Tk License (BSD-style) | https://www.tcl.tk/ |
| Tk | 8.6.15 | TCL | https://www.tcl.tk/ |
| Zstandard |  | BSD 3-Clause OR GPL-2.0 (taken under BSD 3-Clause) | https://facebook.github.io/zstd/ |
| libcurl | 8.18.0 | curl | https://curl.se/ |
| libjpeg-turbo | 3.1.2 | IJG AND BSD-3-Clause AND Zlib | https://libjpeg-turbo.org/ |
| libpng | 1.6.55 | zlib-acknowledgement | http://www.libpng.org/ |
| libtiff | 4.7.0 | HPND | http://www.libtiff.org/ |
| libwebp |  | BSD 3-Clause | https://developers.google.com/speed/webp |
| zlib | 1.3.1 | Zlib | https://zlib.net/ |
| cftime | 1.6.5 | MIT | https://unidata.github.io/cftime/ |
| customtkinter | 5.2.2 | Creative Commons Zero v1.0 Universal | https://customtkinter.tomschimansky.com/ |
| h5py | 3.15.1 | BSD 3-Clause | https://www.h5py.org/ |
| netcdf4 | 1.7.4 | MIT | https://unidata.github.io/netcdf4-python/ |
| numpy | 2.4.2 | BSD-3-Clause | https://numpy.org/ |
| osgeo |  | MIT |  |
| pyhdf | 0.11.6 | MIT | https://github.com/fhs/pyhdf |
| pyproj | 3.7.2 | MIT | https://pyproj4.github.io/pyproj/ |
| rasterio | 1.4.4 | BSD | https://rasterio.readthedocs.io/ |
| scipy | 1.17.1 | BSD-3-Clause | https://scipy.org/ |
| setuptools | 80.10.2 | MIT |  |

## Other files in the build

Present in the build output but not separately attributable components. Most are
Python extension modules belonging to packages already listed above, platform
redistributables, or transitive dependencies of GDAL.

- Windows CRT / API sets (Microsoft redistributable): 42 files
- AWS SDK (transitive, via GDAL cloud drivers): 16 files
- Microsoft Visual C++ runtime: 10 files
- Other libraries and extension modules: 87 files

  Notable DLLs among these:
  - `LIBPQ.dll`
  - `Lerc.dll`
  - `abseil_dll.dll`
  - `aec.dll`
  - `archive.dll`
  - `aws-crt-cpp.dll`
  - `azure-core.dll`
  - `azure-identity.dll`
  - `azure-storage-blobs.dll`
  - `azure-storage-common.dll`
  - `blosc.dll`
  - `bzip2.dll`
  - `cfitsio.dll`
  - `charset.dll`
  - `comerr64.dll`
  - `crc32c.dll`
  - `deflate.dll`
  - `ffi-7.dll`
  - `ffi-8.dll`
  - `ffi.dll`

## Components with obligations beyond attribution

**GEOS** (LGPL-2.1-only)

Distributed as a dynamically linked shared library (DLL) which the user is
free to replace with their own build. This is the form of use the LGPL permits for
software under a different license, and no part of RSDTK is derived from its
source. The library is unmodified from the conda-forge distribution.

**Zstandard** (BSD 3-Clause OR GPL-2.0 (taken under BSD 3-Clause))

Dual-licensed. Used here under the permissive option noted above, which
imposes attribution requirements only. The copyleft option is not exercised.


## Explicitly bundled by the build script

From `build_rsdtk.bat`:

- `--collect-all` rasterio
- `--add-data` %GDAL_DATA_DIR%
- `--add-data` %GDAL_PLUGINS_DIR%
- `--add-data` %PROJ_DATA_DIR%
- `--add-data` docs\RSDTK_User_Guide.html
- `--add-data` ASTER_Converter.py
- `--add-data` Landsat_L1_Converter.py
- `--add-data` Landsat_L2_Converter.py
- `--add-data` Sentinel2_L2A_Converter.py
- `--add-data` Sentinel3_Converter.py
- `--add-data` ECOSTRESS_L1TC_Converter.py
- `--add-data` ECOSTRESS_L2_Converter.py
- `--add-data` VIIRS_Converter.py
- `--add-data` MODIS_Converter.py
- `--add-data` EMIT_L2A_Converter.py
- `--add-data` AIRS_L1B_Converter.py
- `--add-data` MASTER_L1B_Converter.py
- `--add-data` MASTER_L2_Converter.py
- `--add-data` HyTES_L1_Converter.py
- `--add-data` HyTES_L2_Converter.py
- `--add-data` AVIRIS3_L1B_Converter.py
- `--add-data` AVIRIS3_L2A_Converter.py
- `--add-data` AVIRIS_Converter.py
- `--add-data` spectral_tools.py
- `--add-data` subset_tool.py
- `--hidden-import` ASTER_Converter
- `--hidden-import` Landsat_L1_Converter
- `--hidden-import` Landsat_L2_Converter
- `--hidden-import` Sentinel2_L2A_Converter
- `--hidden-import` Sentinel3_Converter
- `--hidden-import` ECOSTRESS_L1TC_Converter
- `--hidden-import` ECOSTRESS_L2_Converter
- `--hidden-import` VIIRS_Converter
- `--hidden-import` MODIS_Converter
- `--hidden-import` EMIT_L2A_Converter
- `--hidden-import` AIRS_L1B_Converter
- `--hidden-import` MASTER_L1B_Converter
- `--hidden-import` MASTER_L2_Converter
- `--hidden-import` HyTES_L1_Converter
- `--hidden-import` HyTES_L2_Converter
- `--hidden-import` AVIRIS3_L1B_Converter
- `--hidden-import` AVIRIS3_L2A_Converter
- `--hidden-import` AVIRIS_Converter
- `--hidden-import` spectral_tools
- `--hidden-import` subset_tool
- `--hidden-import` pyhdf.SD
- `--hidden-import` pyhdf.HDF
- `--hidden-import` pyhdf.V
- `--hidden-import` pyhdf.VS
- `--hidden-import` netCDF4
- `--hidden-import` cftime
- `--hidden-import` h5py
- `--hidden-import` scipy.spatial
- `--hidden-import` scipy.spatial.ckdtree
- `--hidden-import` customtkinter

## Attribution note

MIT and BSD licenses require that the copyright and permission notices be
included with redistributions. Listing the components as above is common
practice for academic software, but the strictly correct approach is to ship the
full license text. The simplest way to do that is to collect each `LICENSE` file
into a `licenses/` folder and have Inno Setup install it alongside the executable.
