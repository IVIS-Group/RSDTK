@echo off
REM ============================================================================
REM  build_rsdtk.bat  —  PyInstaller build script for RSDTK
REM ============================================================================
REM  Usage:
REM      conda activate aster_toolkit
REM      REM      cd /d <path to your RSDTK folder>
REM      build_rsdtk.bat
REM
REM  Output:  dist\RSDTK\RSDTK.exe   (onedir build)
REM  Then:    open RSDTK_Installer.iss in Inno Setup and compile
REM
REM  A command-line build is used rather than a .spec file: PyInstaller 6.x
REM  changed its spec API and the old spec file no longer worked reliably.
REM ============================================================================

echo.
echo ============================================================
echo   Building RSDTK
echo ============================================================
echo.

REM --- Check the conda environment is active -------------------------------
if "%CONDA_PREFIX%"=="" (
    echo ERROR: No conda environment active.
    echo        Run: conda activate aster_toolkit
    exit /b 1
)
echo Environment: %CONDA_PREFIX%
echo.

REM --- Clean previous build ------------------------------------------------
if exist build rmdir /s /q build
if exist dist  rmdir /s /q dist
echo Cleaned build/ and dist/
echo.

REM --- GDAL / PROJ data locations (bundled via runtime_hook.py) -----------
set GDAL_DATA_DIR=%CONDA_PREFIX%\Library\share\gdal
set GDAL_PLUGINS_DIR=%CONDA_PREFIX%\Library\lib\gdalplugins
set PROJ_DATA_DIR=%CONDA_PREFIX%\Library\share\proj
set OPENJP2_DLL=%CONDA_PREFIX%\Library\bin\openjp2.dll

if not exist "%GDAL_DATA_DIR%" echo WARNING: GDAL data not found at %GDAL_DATA_DIR%
if not exist "%GDAL_PLUGINS_DIR%" echo WARNING: gdalplugins not found at %GDAL_PLUGINS_DIR%
if not exist "%PROJ_DATA_DIR%" echo WARNING: PROJ data not found at %PROJ_DATA_DIR%
if not exist "%OPENJP2_DLL%" echo WARNING: openjp2.dll not found at %OPENJP2_DLL%

echo.
echo Running PyInstaller...
echo.

pyinstaller RSDTK.py ^
 --name RSDTK ^
 --onedir ^
 --windowed ^
 --noconfirm ^
 --clean ^
 --runtime-hook runtime_hook.py ^
 --collect-all=rasterio ^
 --collect-submodules=pyproj ^
 --exclude-module PyQt5 ^
 --exclude-module PyQt6 ^
 --exclude-module PySide2 ^
 --exclude-module PySide6 ^
 --exclude-module matplotlib ^
 --add-data "%GDAL_DATA_DIR%;gdal_data" ^
 --add-data "%GDAL_PLUGINS_DIR%;gdalplugins" ^
 --add-data "%PROJ_DATA_DIR%;proj_data" ^
 --add-binary "%OPENJP2_DLL%;." ^
 --add-data "docs\RSDTK_User_Guide.html;." ^
 --hidden-import ASTER_Converter ^
 --hidden-import Landsat_L1_Converter ^
 --hidden-import Landsat_L2_Converter ^
 --hidden-import Sentinel2_L2A_Converter ^
 --hidden-import Sentinel3_Converter ^
 --hidden-import ECOSTRESS_L1TC_Converter ^
 --hidden-import ECOSTRESS_L2_Converter ^
 --hidden-import VIIRS_Converter ^
 --hidden-import MODIS_Converter ^
 --hidden-import EMIT_L2A_Converter ^
 --hidden-import AIRS_L1B_Converter ^
 --hidden-import MASTER_L1B_Converter ^
 --hidden-import MASTER_L2_Converter ^
 --hidden-import HyTES_L1_Converter ^
 --hidden-import HyTES_L2_Converter ^
 --hidden-import AVIRIS3_L1B_Converter ^
 --hidden-import AVIRIS3_L2A_Converter ^
 --hidden-import AVIRIS_Converter ^
 --hidden-import spectral_tools ^
 --hidden-import subset_tool ^
 --hidden-import pyhdf.SD ^
 --hidden-import pyhdf.HDF ^
 --hidden-import pyhdf.V ^
 --hidden-import pyhdf.VS ^
 --hidden-import netCDF4 ^
 --hidden-import cftime ^
 --hidden-import h5py ^
 --hidden-import scipy.spatial ^
 --hidden-import scipy.spatial.ckdtree ^
 --hidden-import customtkinter ^
 --add-data "ASTER_Converter.py;." ^
 --add-data "Landsat_L1_Converter.py;." ^
 --add-data "Landsat_L2_Converter.py;." ^
 --add-data "Sentinel2_L2A_Converter.py;." ^
 --add-data "Sentinel3_Converter.py;." ^
 --add-data "ECOSTRESS_L1TC_Converter.py;." ^
 --add-data "ECOSTRESS_L2_Converter.py;." ^
 --add-data "VIIRS_Converter.py;." ^
 --add-data "MODIS_Converter.py;." ^
 --add-data "EMIT_L2A_Converter.py;." ^
 --add-data "AIRS_L1B_Converter.py;." ^
 --add-data "MASTER_L1B_Converter.py;." ^
 --add-data "MASTER_L2_Converter.py;." ^
 --add-data "HyTES_L1_Converter.py;." ^
 --add-data "HyTES_L2_Converter.py;." ^
 --add-data "AVIRIS3_L1B_Converter.py;." ^
 --add-data "AVIRIS3_L2A_Converter.py;." ^
 --add-data "AVIRIS_Converter.py;." ^
 --add-data "spectral_tools.py;." ^
 --add-data "subset_tool.py;."

if errorlevel 1 (
    echo.
    echo ============================================================
    echo   BUILD FAILED
    echo ============================================================
    exit /b 1
)

echo.
echo ============================================================
echo   BUILD COMPLETE
echo ============================================================
echo.
echo Executable: dist\RSDTK\RSDTK.exe
echo.
echo Next steps:
echo   1. Run dist\RSDTK\RSDTK.exe and confirm it launches
echo   2. Test one granule from a couple of sensors
echo   3. Open RSDTK_Installer.iss in Inno Setup and compile
echo.
