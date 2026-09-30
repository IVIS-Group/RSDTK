# Third-Party Licenses

RSDTK itself is released under the MIT License (see [LICENSE](LICENSE)).

The Windows installer bundles the components below. Running RSDTK from source
installs these through conda or pip instead, in which case their licenses apply
as normal to those installations.

> **TODO before release:** verify this list against the actual build output.
> Run a build, then:
>
>     conda activate aster_toolkit
>     python scripts/collect_licenses.py
>
> That cross-references `dist/RSDTK/_internal/` against the environment's
> package metadata and overwrites this file with verified versions and
> licenses. The list below is the expected set, not a verified one.

| Component | License | Project |
|---|---|---|
| GDAL | MIT | https://gdal.org/ |
| PROJ | MIT | https://proj.org/ |
| rasterio | BSD 3-Clause | https://github.com/rasterio/rasterio |
| NumPy | BSD 3-Clause | https://numpy.org/ |
| SciPy | BSD 3-Clause | https://scipy.org/ |
| pyproj | MIT | https://github.com/pyproj4/pyproj |
| h5py | BSD 3-Clause | https://www.h5py.org/ |
| HDF5 | BSD-style | https://www.hdfgroup.org/ |
| netCDF4-python | MIT | https://github.com/Unidata/netcdf4-python |
| NetCDF-C | BSD 3-Clause | https://www.unidata.ucar.edu/software/netcdf/ |
| pyhdf | MIT | https://github.com/fhs/pyhdf |
| HDF4 | BSD-style | https://www.hdfgroup.org/ |
| CustomTkinter | MIT | https://github.com/TomSchimansky/CustomTkinter |
| OpenJPEG | BSD 2-Clause | https://www.openjpeg.org/ |
| libtiff | libtiff (BSD-style) | http://www.libtiff.org/ |
| Python | PSF License | https://www.python.org/ |

All of the above are permissive licenses compatible with MIT distribution.
None impose copyleft obligations on RSDTK.

## Note on OpenJPEG

The installer bundles `openjp2.dll` and the GDAL plugin directory explicitly, to
enable Sentinel-2 JP2 reading. This is required because PyInstaller does not
detect the GDAL plugin dependency automatically.
