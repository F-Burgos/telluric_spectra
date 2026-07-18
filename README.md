# The Telluric Dataset Project

This repository builds a HARPS atmospheric-spectrum dataset in two stages:

1. Scan selected HARPS science products, extract metadata, and assign a
   practical canonical object name using the existing fast RA-ordered grouping.
2. Stage each object's `e2ds_A` spectra and auxiliary CCF/S1D/BIS products,
   resolve its wavelength/blaze calibrations, and run Naira.

The maintained implementation now lives at the repository root: `phase1/`,
`phase2/`, and `run_pipeline.sh`. Historical copies such as `fase1`, old
`phase2`, `rp1`, and `Results` have been moved into the ignored local archive
`_deprecated/` and are not called by the launcher.

## Environment

The recommended environment manager is `uv`.

```bash
uv sync
```

This creates `.venv/` with the required scientific Python packages. To run the
pipeline through the managed environment:

```bash
uv run ./run_telluric_pipeline.sh \
  --python python
```

Alternatively, activate the environment first:

```bash
source .venv/bin/activate
./run_telluric_pipeline.sh --python python
```

## Data links

The pipeline expects local data entry points under:

```text
Data/science
Data/calib
```

These should usually be symbolic links to the real HARPS data locations. Create
or refresh them with:

```bash
./setup_data_links.sh \
  --science /path/to/HARPS/science \
  --calib /path/to/HARPS/calib
```

Use `--force` to replace existing links:

```bash
./setup_data_links.sh \
  --science /path/to/HARPS/science \
  --calib /path/to/HARPS/calib \
  --force
```

The script only creates symbolic links; it does not copy or move the underlying
data.

## Run

After `Data/science` and `Data/calib` are configured, run from the project root:

```bash
./run_telluric_pipeline.sh
```

By default, generated working outputs are written under:

```text
Output/phase1
Output/spectra
Output/tables
```

The science tree is expected to contain one directory per observing night. By
default Stage 1 opens only `e2ds_A`, CCF-A, S1D-A, and BIS-A FITS products. Use
`--products all` only when a complete FITS inventory is required.

Use `--star NAME` for a single-object smoke run and `--no-fresh` to reuse
existing Stage 1 Parquet files.

Use `--output /path/to/spectra` only when you want Phase 2 somewhere other than
the launch directory's `spectra/` folder.

Telluric PNG previews are disabled by default. To write previews for specific
orders, pass `--plot-orders`:

```bash
./run_telluric_pipeline.sh \
  --plot-orders 63
```

Accepted order selectors are a single order (`63`), a comma-separated list
(`10,20,63`), a range (`50-55`), or `all`.

## Stage 3 dataset build

Stage 2 writes object-oriented cubes. Stage 3 repackages them into a
date-oriented dataset under:

```text
Data/telluric/<NIGHT>/<OBJECT>_<NIGHT>_telluric_cube.fits
```

Run Stage 3 with:

```bash
uv run ./run_stage3_build_dataset.sh --overwrite
```

This also writes:

```text
Output/stage3/telluric_cube_index.parquet
Output/stage3/telluric_cube_index.csv
Output/stage3/telluric_stage3_manifest.json
```

The index table has one row per exposure. Each row points to its Stage 3 cube
and exposure index inside that cube, plus selected atmospheric/header metadata
from the source `e2ds_A` FITS file.

## Principal outputs

- `metadata_raw.parquet`: selected product inventory; `OBJ_ID` preserves the
  original FITS object name and `REL_PATH` is portable relative to the science
  root.
- `metadata.parquet`: filtered working metadata.
- `metadata_final.parquet`: grouped Stage 1 table consumed by Stage 2.
- `Output/spectra/<OBJECT>/<OBJECT>/tell_spec/<OBJECT>_telluric_cube.fits`: one
  object-oriented atmospheric transmission cube, shaped
  `(exposure, order, pixel)`.
- `Data/telluric/<NIGHT>/<OBJECT>_<NIGHT>_telluric_cube.fits`: one
  date-oriented atmospheric transmission cube, shaped
  `(exposure, order, pixel)`.
- `Output/stage3/telluric_cube_index.parquet`: exposure-level index for the
  date-oriented telluric cubes.
- `Output/spectra/<OBJECT>/processing_summary.json`: expected, valid, missing,
  and dimensionally invalid outputs for that object.

The original science header from the first exposure is retained in every
telluric cube FITS file. `TELLCUBE`, `NSPECTRA`, `NORDERS`, and `NPIX` describe
the cube output contract.

The canonical Naira template configuration uses `opt_tmpl = no`, following the
project's current reduction specification. Local HARPS `e2ds_A` inputs have 72
orders, so a two-exposure object has telluric cube shape `(2, 72, 4096)`.
