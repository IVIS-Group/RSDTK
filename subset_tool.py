"""
Universal Subset/Tiling Module
Provides bounding box extraction from multiple vector formats and auto-tiling
for use with satellite data converters (ECOSTRESS, ASTER, VIIRS).

Subsetting modes:
  1. Vector file (Shapefile, KML/KMZ, GeoJSON) - extracts bounding box from geometry
  2. Text file - reads bounding box coordinates from a simple text file
  3. Manual entry - lat/lon bounds specified directly
  4. Auto-tiling - splits full extent into regular tiles
  5. Full extent - no subsetting (default)

Dependencies: GDAL/OGR (via osgeo), numpy
"""

import os
import sys
import re
import zipfile
import tempfile
import numpy as np

try:
    from osgeo import ogr, osr
    ogr.UseExceptions()
except ImportError:
    print("ERROR: GDAL/OGR is required for vector file support.")
    print("Install with: pip install GDAL")
    sys.exit(1)


# ============================================================================
# Bounding box class
# ============================================================================

class BoundingBox:
    """Simple bounding box container with validation and utilities."""

    def __init__(self, min_lat, max_lat, min_lon, max_lon):
        self.min_lat = float(min_lat)
        self.max_lat = float(max_lat)
        self.min_lon = float(min_lon)
        self.max_lon = float(max_lon)
        self._validate()

    def _validate(self):
        if self.min_lat >= self.max_lat:
            raise ValueError(f"min_lat ({self.min_lat}) must be less than "
                             f"max_lat ({self.max_lat})")
        if self.min_lon >= self.max_lon:
            raise ValueError(f"min_lon ({self.min_lon}) must be less than "
                             f"max_lon ({self.max_lon})")
        if not (-90 <= self.min_lat <= 90 and -90 <= self.max_lat <= 90):
            raise ValueError(f"Latitude must be between -90 and 90")
        if not (-180 <= self.min_lon <= 180 and -180 <= self.max_lon <= 180):
            raise ValueError(f"Longitude must be between -180 and 180")

    def intersects(self, other):
        """Check if this bounding box intersects another."""
        return not (self.max_lat < other.min_lat or
                    self.min_lat > other.max_lat or
                    self.max_lon < other.min_lon or
                    self.min_lon > other.max_lon)

    def intersection(self, other):
        """Return the intersection of two bounding boxes, or None if they don't overlap."""
        if not self.intersects(other):
            return None
        return BoundingBox(
            max(self.min_lat, other.min_lat),
            min(self.max_lat, other.max_lat),
            max(self.min_lon, other.min_lon),
            min(self.max_lon, other.max_lon)
        )

    @property
    def width_deg(self):
        return self.max_lon - self.min_lon

    @property
    def height_deg(self):
        return self.max_lat - self.min_lat

    def __repr__(self):
        return (f"BoundingBox(lat=[{self.min_lat:.4f}, {self.max_lat:.4f}], "
                f"lon=[{self.min_lon:.4f}, {self.max_lon:.4f}])")


# ============================================================================
# Bounding box extraction from various formats
# ============================================================================

def bbox_from_vector(filepath):
    """
    Extract bounding box from a vector file (Shapefile, GeoJSON, KML).
    Uses GDAL/OGR to read the geometry and compute the extent.
    Automatically reprojects to WGS84 if needed.

    Parameters:
        filepath: path to vector file (.shp, .geojson, .json, .kml)

    Returns:
        BoundingBox object
    """
    ext = os.path.splitext(filepath)[1].lower()

    # Handle KMZ (zipped KML)
    if ext == '.kmz':
        return _bbox_from_kmz(filepath)

    # Determine driver
    driver_map = {
        '.shp': 'ESRI Shapefile',
        '.geojson': 'GeoJSON',
        '.json': 'GeoJSON',
        '.kml': 'KML',
    }

    driver_name = driver_map.get(ext)
    if driver_name is None:
        raise ValueError(f"Unsupported vector format: {ext}\n"
                         f"Supported: .shp, .geojson, .json, .kml, .kmz")

    # Enable KML driver if needed
    if ext == '.kml':
        drv = ogr.GetDriverByName('KML')
        if drv is None:
            # Try LIBKML as fallback
            drv = ogr.GetDriverByName('LIBKML')
            if drv is None:
                raise RuntimeError("KML driver not available in your GDAL installation.")

    ds = ogr.Open(filepath, 0)  # 0 = read-only
    if ds is None:
        raise IOError(f"Could not open vector file: {filepath}")

    try:
        layer = ds.GetLayer(0)
        if layer is None:
            raise ValueError(f"No layers found in: {filepath}")

        # Get layer extent: (min_x, max_x, min_y, max_y)
        extent = layer.GetExtent()
        min_lon, max_lon, min_lat, max_lat = extent

        # Check if reprojection to WGS84 is needed
        srs = layer.GetSpatialRef()
        if srs is not None:
            wgs84 = osr.SpatialReference()
            wgs84.ImportFromEPSG(4326)
            wgs84.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

            if not srs.IsSame(wgs84):
                srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
                coord_transform = osr.CoordinateTransformation(srs, wgs84)

                # Transform all four corners and get new extent
                corners = [
                    (min_lon, min_lat),
                    (min_lon, max_lat),
                    (max_lon, min_lat),
                    (max_lon, max_lat),
                ]
                transformed_lons = []
                transformed_lats = []
                for x, y in corners:
                    tx, ty, _ = coord_transform.TransformPoint(x, y)
                    transformed_lons.append(tx)
                    transformed_lats.append(ty)

                min_lon = min(transformed_lons)
                max_lon = max(transformed_lons)
                min_lat = min(transformed_lats)
                max_lat = max(transformed_lats)

        bbox = BoundingBox(min_lat, max_lat, min_lon, max_lon)
        print(f"  Vector extent: {bbox}")
        return bbox

    finally:
        ds = None  # Close dataset


def _bbox_from_kmz(filepath):
    """Extract bounding box from a KMZ file (zipped KML)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        with zipfile.ZipFile(filepath, 'r') as z:
            kml_files = [f for f in z.namelist() if f.lower().endswith('.kml')]
            if not kml_files:
                raise ValueError(f"No KML file found inside KMZ: {filepath}")
            z.extract(kml_files[0], tmpdir)
            kml_path = os.path.join(tmpdir, kml_files[0])
            return bbox_from_vector(kml_path)


def bbox_from_textfile(filepath):
    """
    Read bounding box coordinates from a text file.

    Supported formats:
      Format 1 (key=value):
        min_lat=-15.5
        max_lat=-14.0
        min_lon=165.0
        max_lon=167.5

      Format 2 (CSV single line):
        min_lat, max_lat, min_lon, max_lon
        -15.5, -14.0, 165.0, 167.5

      Format 3 (labeled lines):
        min_lat: -15.5
        max_lat: -14.0
        min_lon: 165.0
        max_lon: 167.5

    Parameters:
        filepath: path to text file

    Returns:
        BoundingBox object
    """
    with open(filepath, 'r') as f:
        content = f.read().strip()

    # Try key=value format
    kv_pattern = re.compile(
        r'min_lat\s*[=:]\s*([-\d.]+).*?'
        r'max_lat\s*[=:]\s*([-\d.]+).*?'
        r'min_lon\s*[=:]\s*([-\d.]+).*?'
        r'max_lon\s*[=:]\s*([-\d.]+)',
        re.DOTALL | re.IGNORECASE
    )
    match = kv_pattern.search(content)
    if match:
        min_lat, max_lat, min_lon, max_lon = [float(x) for x in match.groups()]
        bbox = BoundingBox(min_lat, max_lat, min_lon, max_lon)
        print(f"  Text file extent: {bbox}")
        return bbox

    # Try CSV format - look for a line with exactly 4 numbers
    for line in content.split('\n'):
        line = line.strip()
        if not line or line.startswith('#') or line.startswith('min'):
            continue
        parts = re.split(r'[,\s\t]+', line)
        try:
            values = [float(p) for p in parts if p]
            if len(values) == 4:
                min_lat, max_lat, min_lon, max_lon = values
                bbox = BoundingBox(min_lat, max_lat, min_lon, max_lon)
                print(f"  Text file extent: {bbox}")
                return bbox
        except ValueError:
            continue

    raise ValueError(
        f"Could not parse bounding box from: {filepath}\n"
        f"Expected format:\n"
        f"  min_lat = -15.5\n"
        f"  max_lat = -14.0\n"
        f"  min_lon = 165.0\n"
        f"  max_lon = 167.5\n"
        f"Or CSV: -15.5, -14.0, 165.0, 167.5"
    )


def bbox_from_manual(min_lat, max_lat, min_lon, max_lon):
    """
    Create bounding box from manually specified coordinates.

    Parameters:
        min_lat, max_lat, min_lon, max_lon: float values

    Returns:
        BoundingBox object
    """
    bbox = BoundingBox(min_lat, max_lat, min_lon, max_lon)
    print(f"  Manual extent: {bbox}")
    return bbox


def bbox_from_file(filepath):
    """
    Auto-detect file format and extract bounding box.
    Routes to the appropriate reader based on file extension.

    Parameters:
        filepath: path to any supported file format

    Returns:
        BoundingBox object
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")

    ext = os.path.splitext(filepath)[1].lower()

    if ext in ('.shp', '.geojson', '.json', '.kml', '.kmz'):
        return bbox_from_vector(filepath)
    elif ext in ('.txt', '.csv'):
        return bbox_from_textfile(filepath)
    else:
        raise ValueError(
            f"Unsupported file format: {ext}\n"
            f"Supported formats: .shp, .geojson, .json, .kml, .kmz, .txt, .csv"
        )


# ============================================================================
# Auto-tiling
# ============================================================================

class TileGrid:
    """
    Divides a bounding box into regular tiles of specified pixel dimensions.
    """

    def __init__(self, full_bbox, tile_size_pixels=2048, resolution_x_deg=None,
                 resolution_y_deg=None, resolution_m=70.0):
        """
        Parameters:
            full_bbox: BoundingBox of the full extent
            tile_size_pixels: tile dimension in pixels (tiles are square)
            resolution_x_deg: pixel size in degrees (longitude)
            resolution_y_deg: pixel size in degrees (latitude)
            resolution_m: pixel size in meters (used if deg not provided)
        """
        self.full_bbox = full_bbox
        self.tile_size = tile_size_pixels

        # Compute resolution in degrees if not provided
        if resolution_x_deg is None or resolution_y_deg is None:
            center_lat = (full_bbox.min_lat + full_bbox.max_lat) / 2.0
            deg_per_meter = 1.0 / 111320.0
            self.res_y = resolution_m * deg_per_meter
            self.res_x = resolution_m * deg_per_meter / np.cos(np.radians(center_lat))
        else:
            self.res_x = resolution_x_deg
            self.res_y = resolution_y_deg

        # Tile extent in degrees
        self.tile_width_deg = self.tile_size * self.res_x
        self.tile_height_deg = self.tile_size * self.res_y

        # Number of tiles in each direction
        self.n_cols = int(np.ceil(full_bbox.width_deg / self.tile_width_deg))
        self.n_rows = int(np.ceil(full_bbox.height_deg / self.tile_height_deg))
        self.total_tiles = self.n_cols * self.n_rows

    def get_tile_bbox(self, row, col):
        """
        Get the bounding box for a specific tile.

        Parameters:
            row: tile row (0-indexed, from top/north)
            col: tile column (0-indexed, from left/west)

        Returns:
            BoundingBox for this tile
        """
        tile_min_lon = self.full_bbox.min_lon + col * self.tile_width_deg
        tile_max_lon = min(tile_min_lon + self.tile_width_deg, self.full_bbox.max_lon)
        tile_max_lat = self.full_bbox.max_lat - row * self.tile_height_deg
        tile_min_lat = max(tile_max_lat - self.tile_height_deg, self.full_bbox.min_lat)

        return BoundingBox(tile_min_lat, tile_max_lat, tile_min_lon, tile_max_lon)

    def get_all_tiles(self):
        """
        Generator yielding (row, col, tile_bbox) for all tiles.
        """
        for row in range(self.n_rows):
            for col in range(self.n_cols):
                yield row, col, self.get_tile_bbox(row, col)

    def __repr__(self):
        return (f"TileGrid({self.n_cols}x{self.n_rows} tiles, "
                f"{self.tile_size}px, "
                f"tile_size={self.tile_width_deg:.4f}x{self.tile_height_deg:.4f} deg)")


def compute_tile_grid(lat, lon, tile_size_pixels=2048, resolution_m=70.0):
    """
    Compute a tile grid from swath lat/lon arrays.

    Parameters:
        lat, lon: 2D swath coordinate arrays
        tile_size_pixels: tile dimension in pixels
        resolution_m: native pixel resolution in meters

    Returns:
        TileGrid object
    """
    valid = np.isfinite(lat) & np.isfinite(lon) & (lat != 0) & (lon != 0)
    lat_valid = lat[valid]
    lon_valid = lon[valid]

    full_bbox = BoundingBox(
        lat_valid.min(), lat_valid.max(),
        lon_valid.min(), lon_valid.max()
    )

    grid = TileGrid(full_bbox, tile_size_pixels=tile_size_pixels,
                    resolution_m=resolution_m)

    print(f"  Full extent: {full_bbox}")
    print(f"  Tile grid: {grid}")
    print(f"  Total tiles: {grid.total_tiles}")

    return grid


# ============================================================================
# Swath extent helper
# ============================================================================

def get_swath_bbox(lat, lon):
    """
    Compute bounding box from swath lat/lon arrays.

    Parameters:
        lat, lon: 2D arrays

    Returns:
        BoundingBox object
    """
    valid = np.isfinite(lat) & np.isfinite(lon) & (lat != 0) & (lon != 0)
    if not valid.any():
        raise ValueError("No valid geolocation data found.")

    return BoundingBox(
        lat[valid].min(), lat[valid].max(),
        lon[valid].min(), lon[valid].max()
    )


# ============================================================================
# Command-line parsing helpers
# ============================================================================

def parse_subset_args(args):
    """
    Parse subsetting-related command line arguments.

    Recognized arguments:
        --subset <file>             Subset using vector/text file
        --bbox <min_lat> <max_lat> <min_lon> <max_lon>   Manual bounding box
        --tile [size]               Auto-tile with optional pixel size (default 2048)

    Parameters:
        args: list of command line arguments

    Returns:
        dict with keys:
            'mode': 'full' | 'subset' | 'tile'
            'bbox': BoundingBox or None
            'tile_size': int or None
            'subset_file': str or None
    """
    result = {
        'mode': 'full',
        'bbox': None,
        'tile_size': None,
        'subset_file': None,
    }

    i = 0
    while i < len(args):
        arg = args[i].lower()

        if arg == '--subset':
            if i + 1 >= len(args):
                raise ValueError("--subset requires a file path argument")
            result['mode'] = 'subset'
            result['subset_file'] = args[i + 1]
            result['bbox'] = bbox_from_file(args[i + 1])
            i += 2

        elif arg == '--bbox':
            if i + 4 >= len(args):
                raise ValueError("--bbox requires 4 arguments: min_lat max_lat min_lon max_lon")
            try:
                min_lat = float(args[i + 1])
                max_lat = float(args[i + 2])
                min_lon = float(args[i + 3])
                max_lon = float(args[i + 4])
            except ValueError:
                raise ValueError("--bbox arguments must be numeric")
            result['mode'] = 'subset'
            result['bbox'] = bbox_from_manual(min_lat, max_lat, min_lon, max_lon)
            i += 5

        elif arg == '--tile':
            result['mode'] = 'tile'
            # Check if next argument is a number (tile size)
            if i + 1 < len(args) and args[i + 1].isdigit():
                result['tile_size'] = int(args[i + 1])
                i += 2
            else:
                result['tile_size'] = 2048
                i += 1

        else:
            i += 1

    return result


# ============================================================================
# Main - standalone testing
# ============================================================================

if __name__ == '__main__':
    print("Universal Subset/Tiling Module")
    print("=" * 40)

    if len(sys.argv) < 2:
        print(f"\nUsage: python {os.path.basename(__file__)} <vector_file>")
        print(f"\nSupported formats:")
        print(f"  Shapefile:  .shp")
        print(f"  GeoJSON:    .geojson, .json")
        print(f"  KML/KMZ:    .kml, .kmz")
        print(f"  Text/CSV:   .txt, .csv")
        print(f"\nText file format:")
        print(f"  min_lat = -15.5")
        print(f"  max_lat = -14.0")
        print(f"  min_lon = 165.0")
        print(f"  max_lon = 167.5")
        print(f"\nOr as CSV: -15.5, -14.0, 165.0, 167.5")
        sys.exit(0)

    filepath = sys.argv[1]
    bbox = bbox_from_file(filepath)
    print(f"\nExtracted bounding box: {bbox}")
    print(f"  Width:  {bbox.width_deg:.4f} degrees")
    print(f"  Height: {bbox.height_deg:.4f} degrees")

