"""
Generate THIRD_PARTY_LICENSES.md for RSDTK.

Cross-references three sources to work out what the installer actually ships
and under what license:

  1. build_rsdtk.bat        — what the build was told to bundle
  2. dist/RSDTK/_internal/  — what PyInstaller actually bundled (authoritative)
  3. conda-meta/*.json      — the license of each installed package

Run from the repository root, with the build environment active and after a
build has completed:

    conda activate aster_toolkit
    python scripts/collect_licenses.py

Writes THIRD_PARTY_LICENSES.md and prints a summary. Anything it cannot
classify is listed under "Needs manual check" rather than guessed at.
"""

import json
import os
import re
import sys
from collections import OrderedDict

# ---------------------------------------------------------------------------

DIST_DIR = os.path.join("dist", "RSDTK", "_internal")
DIST_DIR_FALLBACK = os.path.join("dist", "RSDTK")
BUILD_BAT = "build_rsdtk.bat"
OUTPUT = "THIRD_PARTY_LICENSES.md"

# Map distributed filenames to the project that owns them. PyInstaller renames
# and flattens a lot, so DLL names rarely match package names.
DLL_TO_PROJECT = {
    "gdal": ("GDAL", "https://gdal.org/"),
    "proj": ("PROJ", "https://proj.org/"),
    "geos": ("GEOS", "https://libgeos.org/"),
    "openjp2": ("OpenJPEG", "https://www.openjpeg.org/"),
    "tiff": ("libtiff", "http://www.libtiff.org/"),
    "libtiff": ("libtiff", "http://www.libtiff.org/"),
    "hdf5": ("HDF5", "https://www.hdfgroup.org/"),
    "hdf": ("HDF4", "https://www.hdfgroup.org/"),
    "mfhdf": ("HDF4", "https://www.hdfgroup.org/"),
    "netcdf": ("NetCDF-C", "https://www.unidata.ucar.edu/software/netcdf/"),
    "sqlite": ("SQLite", "https://www.sqlite.org/"),
    "curl": ("libcurl", "https://curl.se/"),
    "expat": ("Expat", "https://libexpat.github.io/"),
    "zlib": ("zlib", "https://zlib.net/"),
    "zstd": ("Zstandard", "https://facebook.github.io/zstd/"),
    "lz4": ("LZ4", "https://lz4.github.io/lz4/"),
    "jpeg": ("libjpeg-turbo", "https://libjpeg-turbo.org/"),
    "png": ("libpng", "http://www.libpng.org/"),
    "webp": ("libwebp", "https://developers.google.com/speed/webp"),
    "xml2": ("libxml2", "https://gitlab.gnome.org/GNOME/libxml2"),
    "openssl": ("OpenSSL", "https://www.openssl.org/"),
    "crypto": ("OpenSSL", "https://www.openssl.org/"),
    "ssl": ("OpenSSL", "https://www.openssl.org/"),
    "tcl": ("Tcl", "https://www.tcl.tk/"),
    "tk": ("Tk", "https://www.tcl.tk/"),
    "python3": ("Python", "https://www.python.org/"),
    "libblas": ("OpenBLAS", "https://www.openblas.net/"),
    "openblas": ("OpenBLAS", "https://www.openblas.net/"),
    "lapack": ("LAPACK", "https://www.netlib.org/lapack/"),
}

# Python packages we expect to see, with their project URLs
PY_PACKAGE_URLS = {
    "numpy": "https://numpy.org/",
    "scipy": "https://scipy.org/",
    "rasterio": "https://rasterio.readthedocs.io/",
    "pyproj": "https://pyproj4.github.io/pyproj/",
    "h5py": "https://www.h5py.org/",
    "netcdf4": "https://unidata.github.io/netcdf4-python/",
    "pyhdf": "https://github.com/fhs/pyhdf",
    "customtkinter": "https://customtkinter.tomschimansky.com/",
    "darkdetect": "https://github.com/albertosottile/darkdetect",
    "packaging": "https://github.com/pypa/packaging",
    "certifi": "https://github.com/certifi/python-certifi",
    "affine": "https://github.com/rasterio/affine",
    "attrs": "https://www.attrs.org/",
    "click": "https://palletsprojects.com/p/click/",
    "cligj": "https://github.com/mapbox/cligj",
    "snuggs": "https://github.com/mapbox/snuggs",
    "pyparsing": "https://github.com/pyparsing/pyparsing",
    "cftime": "https://unidata.github.io/cftime/",
    "pillow": "https://python-pillow.org/",
}


def find_dist_dir():
    for d in (DIST_DIR, DIST_DIR_FALLBACK):
        if os.path.isdir(d):
            return d
    return None


def scan_dist(dist):
    """Return (dll_projects, py_packages) actually present in the build."""
    dll_projects = {}
    py_packages = set()
    unclassified = set()

    for entry in os.listdir(dist):
        path = os.path.join(dist, entry)
        low = entry.lower()

        if low.endswith(".dll") or low.endswith(".pyd"):
            stem = re.sub(r"[-_.]?\d.*$", "", os.path.splitext(low)[0])
            stem = stem.lstrip("lib")
            matched = False
            for key, (proj, url) in DLL_TO_PROJECT.items():
                if stem.startswith(key) or key in stem:
                    dll_projects[proj] = url
                    matched = True
                    break
            if not matched and not low.startswith("_"):
                unclassified.add(entry)

        elif os.path.isdir(path):
            name = low.replace("-", "_")
            if name in ("gdal_data", "gdalplugins", "proj_data", "share",
                        "tcl", "tk", "tcl8", "tk8"):
                continue
            py_packages.add(name)

    return dll_projects, py_packages, unclassified


def read_conda_licenses():
    """Map package name -> (version, license) from conda-meta JSON files."""
    prefix = os.environ.get("CONDA_PREFIX")
    if not prefix:
        return {}
    meta = os.path.join(prefix, "conda-meta")
    if not os.path.isdir(meta):
        return {}

    out = {}
    for fn in os.listdir(meta):
        if not fn.endswith(".json"):
            continue
        try:
            with open(os.path.join(meta, fn), "r", encoding="utf-8") as fh:
                d = json.load(fh)
        except (OSError, json.JSONDecodeError):
            continue
        name = (d.get("name") or "").lower()
        if name:
            out[name] = (d.get("version", "?"), d.get("license") or "")
    return out


def read_pip_licenses():
    """Map package name -> (version, license) from installed dist-info."""
    out = {}
    try:
        from importlib import metadata
    except ImportError:
        return out

    for dist in metadata.distributions():
        try:
            name = (dist.metadata["Name"] or "").lower()
        except Exception:
            continue
        if not name:
            continue
        lic = dist.metadata.get("License") or ""
        if not lic or len(lic) > 80:
            # Some packages dump the whole license text here; prefer classifiers
            classifiers = dist.metadata.get_all("Classifier") or []
            for c in classifiers:
                if c.startswith("License ::"):
                    lic = c.split("::")[-1].strip()
                    break
        out[name] = (dist.version, lic)
    return out


def parse_build_bat():
    """Extract explicitly bundled items from build_rsdtk.bat, if present."""
    if not os.path.exists(BUILD_BAT):
        return []
    with open(BUILD_BAT, "r", encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    items = []
    for m in re.finditer(r"--collect-all[= ]+([\w.]+)", text):
        items.append(("collect-all", m.group(1)))
    for m in re.finditer(r'--add-data[= ]+"?([^";\s]+)', text):
        items.append(("add-data", m.group(1)))
    for m in re.finditer(r"--hidden-import[= ]+([\w.]+)", text):
        items.append(("hidden-import", m.group(1)))
    return items


def lookup_license(name, conda, pip):
    key = name.lower().replace("-", "_")
    for table in (conda, pip):
        if key in table:
            return table[key]
        # try a few normalizations
        for alt in (key.replace("_", "-"), key.replace("4", ""), "lib" + key):
            if alt in table:
                return table[alt]
    return None


def main():
    dist = find_dist_dir()
    conda = read_conda_licenses()
    pip = read_pip_licenses()
    bat_items = parse_build_bat()

    print(f"conda-meta packages found: {len(conda)}")
    print(f"installed Python distributions: {len(pip)}")
    print(f"build_rsdtk.bat entries: {len(bat_items)}")

    if dist:
        print(f"scanning build output: {dist}")
        dll_projects, py_packages, unclassified = scan_dist(dist)
        print(f"  native libraries identified: {len(dll_projects)}")
        print(f"  Python packages bundled: {len(py_packages)}")
        print(f"  unclassified files: {len(unclassified)}")
    else:
        print("WARNING: no dist/ folder found. Run a build first for an")
        print("         authoritative list. Falling back to the expected set.")
        dll_projects = {p: u for p, u in DLL_TO_PROJECT.values()}
        py_packages = set(PY_PACKAGE_URLS)
        unclassified = set()

    # --- assemble rows ---
    rows = []
    needs_check = []

    for proj in sorted(dll_projects):
        url = dll_projects[proj]
        info = lookup_license(proj, conda, pip)
        if info:
            rows.append((proj, info[0], info[1], url))
        else:
            rows.append((proj, "", "", url))
            needs_check.append(proj)

    for pkg in sorted(py_packages):
        info = lookup_license(pkg, conda, pip)
        url = PY_PACKAGE_URLS.get(pkg, "")
        if info:
            rows.append((pkg, info[0], info[1], url))
            if not info[1]:
                needs_check.append(pkg)
        else:
            rows.append((pkg, "", "", url))
            needs_check.append(pkg)

    # de-duplicate by name, preferring rows with a license
    merged = OrderedDict()
    for name, ver, lic, url in rows:
        k = name.lower()
        if k not in merged or (not merged[k][2] and lic):
            merged[k] = (name, ver, lic, url)

    # --- write ---
    with open(OUTPUT, "w", encoding="utf-8") as fh:
        fh.write("# Third-Party Licenses\n\n")
        fh.write("RSDTK is released under the MIT License "
                 "(see [LICENSE](LICENSE)).\n\n")
        fh.write("The Windows installer bundles the components below. Running "
                 "RSDTK from source\ninstalls these via conda or pip instead, "
                 "in which case their licenses apply to\nthose installations "
                 "in the normal way.\n\n")
        if dist:
            fh.write(f"Generated from the build output in `{dist}` and the "
                     "build environment's package\nmetadata, by "
                     "`scripts/collect_licenses.py`.\n\n")
        else:
            fh.write("> **Note:** generated without a build present, so this "
                     "is the expected set rather\n> than a verified one. "
                     "Re-run after a build.\n\n")

        fh.write("| Component | Version | License | Project |\n")
        fh.write("|---|---|---|---|\n")
        for name, ver, lic, url in merged.values():
            lic_cell = lic if lic else "**check**"
            url_cell = url if url else ""
            fh.write(f"| {name} | {ver} | {lic_cell} | {url_cell} |\n")

        if needs_check:
            fh.write("\n## Needs manual check\n\n")
            fh.write("No license could be determined automatically for these. "
                     "Look them up and fill\nthem into the table above.\n\n")
            for n in sorted(set(needs_check)):
                fh.write(f"- {n}\n")

        if unclassified:
            fh.write("\n## Unrecognized files in the build\n\n")
            fh.write("Present in the build output but not matched to a known "
                     "project. Most will be\nPython extension modules "
                     "belonging to packages already listed.\n\n")
            for f in sorted(unclassified)[:60]:
                fh.write(f"- `{f}`\n")
            if len(unclassified) > 60:
                fh.write(f"- ... and {len(unclassified) - 60} more\n")

        if bat_items:
            fh.write("\n## Explicitly bundled by the build script\n\n")
            fh.write("From `build_rsdtk.bat`:\n\n")
            for kind, val in bat_items:
                fh.write(f"- `--{kind}` {val}\n")

        fh.write("\n## Attribution note\n\n")
        fh.write("MIT and BSD licenses require that the copyright and "
                 "permission notices be\nincluded with redistributions. "
                 "Listing the components as above is common\npractice for "
                 "academic software, but the strictly correct approach is to "
                 "ship the\nfull license text. The simplest way to do that is "
                 "to collect each `LICENSE` file\ninto a `licenses/` folder "
                 "and have Inno Setup install it alongside the "
                 "executable.\n")

    print(f"\nWrote {OUTPUT}")
    print(f"  {len(merged)} components listed")
    if needs_check:
        print(f"  {len(set(needs_check))} need a manual license lookup")
    return 0


if __name__ == "__main__":
    sys.exit(main())
