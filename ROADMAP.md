# RSDTK Roadmap

Planning document for the Remote Sensing Data Toolkit: dissemination strategy and
staged development work.

**Current version:** v1.2
**Status:** Preparing for public release

---

## Contents

- [Dissemination plan](#dissemination-plan)
- [Timeline](#timeline)
- [Software: Gate 1 — before the repository goes public](#software-gate-1--before-the-repository-goes-public)
- [Software: Gate 2 — before paper submission](#software-gate-2--before-paper-submission)
- [Software: Gate 3 — post-release and community-driven](#software-gate-3--post-release-and-community-driven)
- [Dataset gaps](#dataset-gaps)
- [Notes on venue choices](#notes-on-venue-choices)

---

## Dissemination plan

Three parallel efforts, each pointing at the same public release.

### Prong 1 — Open release

Publish RSDTK on the lab GitHub organization with a tagged release and a Zenodo
DOI wired to the GitHub release integration, so every version is archived and
independently citable. A concept DOI always resolves to the latest version.

This is the foundation the other two prongs reference.

### Prong 2 — Software description paper

A paper in which the toolkit itself is the subject, not the instrument.

Candidate venues, to be settled once the University of Pittsburgh APC
renegotiation concludes:

| Venue | Publisher | Notes |
|---|---|---|
| *Environmental Modelling & Software* | Elsevier | Software central to scope; potentially APC-covered |
| *Earth and Space Science* | AGU / Wiley | Methods, instruments, data and software; potentially APC-covered |
| *Journal of Open Research Software* | Ubiquity | Has a dedicated Software Metapaper article type |
| *Geoscientific Model Development* | Copernicus | Excellent fit, but not APC-covered |

Deliberately **not** *Computers & Geosciences*: its article types require either an
original research contribution or a third-party software review, neither of which
fits a self-authored tool description. C&G remains a strong target for a later
research paper that uses RSDTK as apparatus.

Deliberately **not** *Journal of Open Source Software*: the six-month public
development history requirement and the automated-test and public-code-review
obligations carry more overhead than the format justifies given competing
commitments. This decision can be revisited later at no cost, since the
repository will accumulate exactly the history JOSS looks for.

### Prong 3 — Teaching dissemination

**SERC / Teach the Earth activities.** Submit one or two teaching activities to
the Teach the Earth portal. Activities are accepted at any time and reviewed on a
rolling basis; those passing review are tagged and rank higher in portal search
results. Reviewers score five elements: scientific accuracy, alignment of goals /
activity / assessment, pedagogic effectiveness, robustness, and completeness of
the activity description.

Candidate activities, in order of how well they exploit RSDTK's distinctive
capability:

1. **Multi-sensor thermal comparison of an active volcano.** The same target
   through ASTER, ECOSTRESS, Landsat and MODIS. Learning objective: sensor
   selection as a scientific decision, and the resolution / revisit trade-off.
   Infeasible to assign without a harmonising tool, which is the argument.
2. **Spectral resolution and what it buys.** EMIT hyperspectral versus Landsat
   multispectral over one surface; students resample EMIT to Landsat bands using
   the match-sensor feature and identify which features survive.
3. **The fill-value trap.** Students compute NDVI on data with unhandled fill
   values, discover the spurious values of 1, then diagnose and fix it. Small,
   self-contained, and teaches data scepticism by design rather than by warning.
   Best candidate for a first submission.

**Critical design constraint:** bundle pre-processed output files with every
activity. This decouples adoption of the activity from installation of the
software, lets instructors on any platform evaluate and run it, and makes direct
RSDTK use an optional extension rather than a prerequisite. Reviewers assess
robustness, and an activity requiring Windows-only software scores poorly.

**Education paper.** A Curriculum and Instruction submission to the *Journal of
Geoscience Education*, which requires evidence of effectiveness for new teaching
materials. Taylor & Francis is not among Pitt's APC agreements, but JGE is a
subscription journal so publication carries no charge if open access is not
elected.

Evidence to gather while teaching:

- Pre/post self-efficacy survey on handling multi-sensor data
- Time to first analysis-ready product, ideally against a prior cohort who
  processed data manually
- Number of distinct sensors each student engaged with (the core feasibility
  claim)
- Rubric-scored final projects
- Error audit: frequency of data-preparation mistakes (wrong scale factors, fill
  values treated as valid, CRS mismatches). The most objective available measure,
  and specific to what RSDTK prevents

Framing note: cognitive load theory provides a citable structure, where
extraneous load from file handling crowds out germane load on radiative transfer
and spectral interpretation. A reviewer will ask whether abstracting away
preprocessing deprives students of necessary understanding. Build the answer into
the course design rather than deflecting it: use RSDTK early to reach science
quickly, then have students implement one conversion by hand, and measure whether
that sequence outperforms either extreme.

---

## Timeline

| Period | Work |
|---|---|
| **Immediate** | IRB contact; secrets and path audit; license; README; repository public |
| **October 2026** | Zenodo DOI; tagged v1.2 release; CHANGELOG; CITATION.cff |
| **Autumn term** | Student evidence collection; capture course materials while teaching; direct user reports to GitHub Issues |
| **Nov–Dec 2026** | Cross-platform source install; automated tests; provenance metadata |
| **Winter break** | Software paper draft |
| **Spring 2027** | Paper submission; SERC activity submissions; APC decision resolved |
| **Summer 2027 onward** | Education paper, once assessment evidence is sufficient |

### The one hard deadline

**Contact the Pitt Human Research Protection Office about IRB now.** Classroom
education research is frequently granted exempt status, but the determination
generally needs to predate data collection if the results are to be published,
and students are using RSDTK this term. Even where exemption applies, having the
determination on file before collecting with publication intent matters, and JGE
will ask.

This is the only item in the plan with a deadline that can be permanently missed.

---

## Software: Gate 1 — before the repository goes public

Items that are awkward or impossible to correct after the first push.

### Security and hygiene

- [ ] **Absolute path audit.** Grep the codebase for hard-coded paths; several are
      expected given the development layout under
      `D:\University\...\RSDTK`.
- [ ] **Credential audit.** Check for EarthData credentials, API tokens, or
      `.netrc` references. Anything committed persists in git history even after
      deletion, so this must precede the first push.
- [ ] **`.gitignore`** covering `__pycache__/`, `build/`, `dist/`, `*.spec`,
      `*.pyc`, and test data directories.
- [ ] **Do not commit the installer binary.** A 100+ MB executable permanently
      bloats the repository. Attach it to a GitHub Release instead.

### Licensing

- [ ] **`LICENSE`** — MIT recommended. The dependency stack is permissively
      licensed throughout, so nothing constrains the choice.
- [ ] **`THIRD_PARTY_LICENSES.md`** — the installer redistributes GDAL and related
      binaries, whose licenses should be included.

### Correctness of what is claimed

- [ ] **Strip debugging artifacts.** Any diagnostic log lines, `apply_patches.py`,
      and remaining scratch scripts.
- [ ] **Decide on `ACOLITE_L2W_Converter.py`.** Built as a standalone for ACOLITE
      Level-2 Water NetCDF and never integrated into the GUI, given how niche the
      format is. Either ship it in a `standalone/` directory with a short note, or
      leave it out. Shipping it unexplained will generate questions.
- [ ] **Spot-check the user guide against the code**, sensor and product tables in
      particular. The guide is the source for the README table, so an error there
      propagates.

*Checked and already resolved during development, recorded here so they are not
re-investigated: the AVIRIS-3 resolution field works (resampling added in v1.0);
ASTER V003 including AST_07M is implemented (v1.0); light mode appearance was
polished in v1.0; `_organize_output_files()` was removed in v1.1 because it moved
original input files, and each converter creates its own `GeoTIFF/` and `ENVI/`
subdirectories internally.*

### Repository files

- [ ] **`README.md`** — the highest-value file in the repository. Should contain:
      what the toolkit does; the sensor and product support table; installation
      for both the Windows installer and from source; a quickstart; a GUI
      screenshot; citation information; license; and a note that automated tests
      are not yet in place.
- [ ] **`requirements.txt`** and **`environment.yml`** with pinned versions
      exported from the `aster_toolkit` environment.
- [ ] **`CHANGELOG.md`** — buildable from the user guide's "What's New" sections,
      with dates.
- [ ] **`CITATION.cff`** — renders the GitHub "Cite this repository" button and
      supplies Zenodo with correct metadata.
- [ ] **`CONTRIBUTING.md`** and a bug-report issue template prompting for sensor,
      product type, platform, and log output.

### Git history

Push genuine commit history if it exists. If development was not under version
control, start clean and document the timeline honestly in the CHANGELOG. Do not
fabricate backdated commits: committer dates are exposed separately from author
dates and are straightforward to inspect.

---

## Software: Gate 2 — before paper submission

### Cross-platform source installation

The single most important item in this tier. A reviewer must be able to run the
software, and a Windows-only tool will be flagged at any software venue. Platform
installers are not required, but a working `pip install -e .` path is, verified
far enough on macOS and Linux that the GUI launches and at least one converter
completes.

This is also the largest adoption barrier for colleagues generally.

### Automated tests

Not coverage for its own sake. Target the failure modes observed repeatedly
during development:

- [ ] Parameter threading through the full call chain — `bbox`, `spectral_opts`
      and `overwrite` reaching the write functions, which has silently broken
      more than once
- [ ] Scale factors applied correctly per product (AST_05/07 × 0.001,
      AST_08 × 0.1, AST_09T per-band UCC)
- [ ] Fill-value handling: NaN in GeoTIFF, −9999 in ENVI, and no propagation of
      zeros as valid data
- [ ] Spectral band masks selecting the bands they claim, including the
      nanometer / micrometer unit conversion
- [ ] Bounding-box clipping producing the expected output dimensions

Use synthetic fixtures generated within the tests rather than shipping real
granules.

### Reproducibility

- [ ] **Provenance metadata in outputs.** RSDTK version, processing timestamp and
      parameters written into GeoTIFF tags and the ENVI header. Inexpensive, and
      reviewers will ask.
- [ ] **Batch configuration files** (YAML or JSON). Makes processing recipes
      reproducible and scriptable, and gives users something specific to cite in
      a methods section.
- [ ] **A runnable demonstration.** Either a small clipped sample committed to the
      repository, or a script fetching a named granule from EarthData. Without
      this a reviewer cannot verify that anything works.

### Highest-priority dataset gaps

- [ ] **MODIS L1B radiance** (MOD021KM / MOD02HKM / MOD02QKM). The largest gap
      relative to the intended audience: bands 21 and 22 at 3.9 µm are the
      standard MODIS thermal anomaly channels.
- [ ] **Sentinel-2 L1C.** TOA reflectance for users running their own atmospheric
      correction, plus B10 (1375 nm cirrus) which exists only in L1C, and the 60 m
      band group (B01, B09) currently unsupported.
- [ ] **VIIRS NOAA-20/21 prefixes** (`VJ102`/`VJ103`, `VJ202`/`VJ203`). Close to
      free and an awkward gap for a user to discover first.

### Documentation

- [ ] Publish the HTML user guide in the repository as a `docs/` folder or a
      GitHub Pages build, not only as a PDF.

---

## Software: Gate 3 — post-release and community-driven

Publishing this tier as a GitHub milestone or roadmap gives the community
something concrete to react to, which generates the issue traffic that
strengthens both papers.

### Functionality

- [ ] Arbitrary EPSG output rather than Geographic / UTM only
- [ ] Bilinear and cubic resampling exposed in the GUI (already a v1.3 candidate)
- [ ] Cloud Optimized GeoTIFF output
- [ ] Multiprocessing, primarily to address VIIRS processing time (~15 min/file)
- [ ] Log-to-file, so users send a file rather than a screenshot
- [ ] Per-file progress rather than an indeterminate progress bar

### Reprojection quality

- [ ] Improve KD-tree reprojection for Sentinel-3 and other swath data: bilinear
      interpolation, UTM output rather than geographic WGS-84, and cos(lat)
      correction for longitude distances

---

## Dataset gaps

Full list, beyond the Gate 2 priorities above. The pattern to watch is a
supported mission whose *other* products fail, which reads as a bug rather than
as scope.

| Sensor | Missing product | Rationale |
|---|---|---|
| MODIS | MOD14 thermal anomalies | Companion to MOD021KM for thermal work |
| VIIRS | VNP14 / VNP14IMG active fire | Successor to MODVOLC-style monitoring |
| VIIRS | VNP46 Black Marble | Nighttime lights, eruption glow |
| ECOSTRESS | Tiled collection (L2T_LSTE, L3T, L4T) | Already on UTM grids; if tiles supersede swath, users will hit "no files found" |
| ASTER | AST_09, AST_09XT | VNIR/SWIR surface radiance; only 09T TIR supported |
| ASTER | ASTER GED | Gridded emissivity, widely used as a TIR prior |
| EMIT | L2B mineralogy, L1B radiance | L2B MIN is the flagship EMIT product |
| Sentinel-3 | SLSTR L2 FRP | Fire radiative power, direct volcanic application |
| Landsat | Collection 2 Level-3 (DSWE, fSCA, burned area) | |
| Landsat | SLC-off gap flagging for L7 post-2003 | Otherwise stripes pass through as valid data |
| AVIRIS-3/5 | L2A surface reflectance | Pending data access via the AVIRIS portal |
| Sounders | CrIS | Operational successor to AIRS, same use case |

Longer-term sensor additions worth considering if community interest supports
them: GOES-R ABI and Himawari AHI for high-temporal volcanic monitoring, MSG
SEVIRI, SDGSAT-1 TIS, PACE OCI, and TROPOMI for volcanic SO₂.

---

## Notes on venue choices

**Why a paper rather than a DOI alone.** A Zenodo DOI is a citation mechanism,
not a publication. It provides no discovery pathway: nobody finds a tool by
browsing Zenodo. A software description paper is indexed, appears in Google
Scholar, and gives the toolkit a citable anchor that does not require first
writing a separate research paper that happens to use it.

**Why the three prongs reinforce each other.** SERC activities are evidence of
use and dissemination, which is what the JGE paper requires and what a software
paper's statement of need rests on. The activities are not a detour from the
papers; they are the substrate for them. Likewise, the GitHub issue history
generated by real users is the community-feedback evidence that strengthens any
software submission.

**On soliciting community suggestions.** Planned for after release. Publishing
Gate 3 as an open roadmap makes that solicitation concrete rather than
open-ended, and gives potential contributors an obvious entry point.
