"""
Runtime hook for PyInstaller — sets GDAL/PROJ environment variables
so that rasterio and pyproj can find their data files in the bundled app.
"""

import os
import sys

if getattr(sys, 'frozen', False):
    bundle_dir = sys._MEIPASS

    # GDAL support data (coordinate system tables, etc.)
    gdal_data = os.path.join(bundle_dir, 'gdal_data')
    if os.path.isdir(gdal_data):
        os.environ['GDAL_DATA'] = gdal_data

    # GDAL format plugins (needed for Sentinel-2 JP2 via JP2OpenJPEG)
    gdal_plugins = os.path.join(bundle_dir, 'gdalplugins')
    if os.path.isdir(gdal_plugins):
        os.environ['GDAL_DRIVER_PATH'] = gdal_plugins

    # PROJ database — location varies between builds, so try each in turn
    for proj_candidate in ('proj_data', 'proj', 'share/proj'):
        proj_path = os.path.join(bundle_dir, proj_candidate)
        if os.path.isdir(proj_path):
            os.environ['PROJ_LIB'] = proj_path
            break
