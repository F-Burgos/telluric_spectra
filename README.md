# The Telluric Dataset Project

This repository builds a HARPS atmospheric-spectrum dataset in two stages:

1. Scan selected HARPS science products, extract metadata, and assign a
   practical canonical object name using the existing fast RA-ordered grouping.
2. Stage each object's `e2ds_A` spectra and auxiliary CCF/S1D/BIS products,
   resolve its wavelength/blaze calibrations, and run Naira.

The maintained implementation is under `p_apoyo/p_apoyo`. Historical root-level
copies such as `fase1`, `phase2`, `rp1`, and `Results` have been moved into the
ignored local archive `_deprecated/` and are not called by the launcher.

## Run

Provide pointers to the science and calibration trees. By default, Stage 2 writes
object folders under `./spectra`, relative to the directory where you launch the
pipeline:

```bash
./run_telluric_pipeline.sh \
  --science /path/to/HARPS/science \
  --calib /path/to/HARPS/calib \
  --work-dir /path/to/work/phase1 \
  --tables-path /path/to/output/tables
```

The science tree is expected to contain one directory per observing night. By
default Stage 1 opens only `e2ds_A`, CCF-A, S1D-A, and BIS-A FITS products. Use
`--products all` only when a complete FITS inventory is required.

Use `--star NAME` for a single-object smoke run and `--no-fresh` to reuse
existing Stage 1 Parquet files.

Use `--output /path/to/spectra` only when you want Phase 2 somewhere other than
the launch directory's `spectra/` folder.

## Principal outputs

- `metadata_raw.parquet`: selected product inventory; `OBJ_ID` preserves the
  original FITS object name and `REL_PATH` is portable relative to the science
  root.
- `metadata.parquet`: filtered working metadata.
- `metadata_final.parquet`: grouped Stage 1 table consumed by Stage 2.
- `spectra/<OBJECT>/<OBJECT>/tell_spec/*_telluric.fits`: one atmospheric
  transmission matrix per accepted exposure, shaped `(order, pixel)`.
- `spectra/<OBJECT>/processing_summary.json`: expected, valid, missing,
  and dimensionally invalid outputs for that object.

The original science header is retained in every telluric FITS file. `TELLMAT`,
`NORDERS`, and `NPIX` describe the matrix output contract.

The canonical Naira template configuration uses `opt_tmpl = no`, following the
project's current reduction specification. Local HARPS `e2ds_A` inputs have 72
orders, so their telluric matrices have shape `(72, 4096)`.
