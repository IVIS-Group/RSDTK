# Contributing to RSDTK

Contributions are welcome, and that includes bug reports. If you hit a problem
with a sensor or product, reporting it is a genuine contribution — several
converters exist in their current form because someone said "this doesn't work
with my data."

## Reporting a bug

Open an issue using the bug report template. The most useful things to include:

- **Sensor and product** — e.g. "ECOSTRESS L2 LSTE", "Landsat 7 Level-2"
- **The full Processing Log output.** Scroll to capture the whole traceback, not
  just the last line. This is by far the most useful single item.
- **How the input folder is arranged** — one subfolder per granule, or flat
- **A filename or two** from the input data. Product detection works from
  filename patterns, so a mismatch is often the cause.
- **Your platform** and whether you used the installer or ran from source

Please do not attach granules. They are large and usually unnecessary. If a file
genuinely seems malformed, the filename and a description of the internal
structure is normally enough to work from.

## Requesting a sensor or product

Open an issue and include:

- Which mission, instrument, and product level
- Where the data is distributed from
- A representative filename
- Roughly what the internal structure looks like, if you know. For HDF or
  NetCDF, the output of `h5dump -n file.h5`, `ncdump -h file.nc`, or a short
  `pyhdf` dataset listing is ideal.
- What you would do with it, which helps decide sensible defaults

Known gaps are already listed in [`ROADMAP.md`](ROADMAP.md). A comment on an
existing item is as useful as a new issue, since priorities are driven partly by
how many people ask.

## Contributing code

### Setting up

```bash
git clone https://github.com/<ORG>/RSDTK.git
cd RSDTK
conda env create -f environment.yml
conda activate rsdtk
python RSDTK.py
```

### Adding a converter

Existing converters follow a consistent shape, and matching it makes review
easier:

- `discover_*()` — find and group input files, returning granule identifiers
- `read_*_data()` — read arrays and metadata, apply scaling and fill values
- `compute_output_grid()` — build the output grid, honoring `bbox` if given
- `process_granule()` — orchestrate one granule
- `process_directory()` — iterate granules, and serve as the GUI entry point
- `write_output()` — write GeoTIFF and/or ENVI

`process_directory()` is what `RSDTK.py` imports, so its signature is the
integration point.

### Conventions worth knowing

These are settled decisions, arrived at the hard way:

- **Fill values.** NaN in GeoTIFF, `-9999` in ENVI. Never zero, which produces
  spurious index values such as NDVI = 1.
- **ENVI headers.** Write the binary and header directly with NumPy rather than
  via rasterio's ENVI driver, which discards wavelength metadata.
- **Swath reprojection.** cKDTree with inverse distance weighting, and
  zero-initialized accumulation rather than NaN-initialized, which otherwise
  propagates NaN through the output.
- **GLT orthorectification** where the product provides a geolocation lookup
  table. Considerably faster than KD-tree and more accurate.
- **NetCDF reads.** Call `set_auto_maskandscale(False)` before reading, or
  scale and offset get applied twice.
- **Parameter threading.** New parameters such as `bbox` or `spectral_opts` must
  be added to every signature in the chain, from `process_directory()` down to
  the write functions. Missing one link fails silently, with the option simply
  having no effect. This has happened more than once; please check the whole
  chain.
- **QA layers** are written as separate files, not appended to science bands.
- **Wavelength units.** Micrometers in ENVI headers for TIR, nanometers for
  VSWIR. `spectral_tools.py` works in nanometers throughout, so convert on the
  way in and out.

### Before opening a pull request

- Confirm the GUI still launches and the sensor you touched still processes
- Run the test suite if one exists by the time you read this
- Check that any new option actually takes effect end to end, rather than being
  accepted and ignored
- Keep to US English spelling in code comments and user-facing text

## Questions

Open an issue with the question label. If something in the documentation was
unclear enough to prompt a question, that is worth knowing about too.
