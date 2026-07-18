"""Build the date-oriented telluric cube dataset.

Stage 2 writes one cube per object.  Stage 3 repackages those cubes by
observing night:

    Data/telluric/<NIGHT>/<OBJECT>_<NIGHT>_telluric_cube.fits

It also writes an exposure-level index table.  Each table row points to the
night cube plus the exposure index inside that cube, and records selected
metadata extracted from the original e2ds FITS header.
"""

from __future__ import annotations

import argparse
import json
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits


HEADER_KEYS = {
    "ESO TEL AIRM END": "ESO_TEL_AIRM_END",
    "ESO TEL AIRM START": "ESO_TEL_AIRM_START",
    "ESO TEL ALT": "ESO_TEL_ALT",
    "ESO TEL AMBI FWHM END": "ESO_TEL_AMBI_FWHM_END",
    "ESO TEL AMBI FWHM START": "ESO_TEL_AMBI_FWHM_START",
    "ESO TEL AMBI PRES END": "ESO_TEL_AMBI_PRES_END",
    "ESO TEL AMBI PRES START": "ESO_TEL_AMBI_PRES_START",
    "ESO TEL AMBI RHUM": "ESO_TEL_AMBI_RHUM",
    "ESO TEL AMBI TEMP": "ESO_TEL_AMBI_TEMP",
    "ESO TEL AMBI WINDDIR": "ESO_TEL_AMBI_WINDDIR",
    "ESO TEL AMBI WINDSP": "ESO_TEL_AMBI_WINDSP",
    "ESO DRS BERV": "ESO_DRS_BERV",
    "ESO DRS SPE EXT SN60": "ESO_DRS_SPE_EXT_SN60",
    "ESO TEL MOON RA": "ESO_TEL_MOON_RA",
    "ESO TEL MOON DEC": "ESO_TEL_MOON_DEC",
    "DATE-OBS": "DATE_OBS",
    "MJD-OBS": "MJD_OBS",
    "RA": "RA",
    "DEC": "DEC",
}


def safe_name(value: object) -> str:
    """Return a filesystem-friendly name while keeping object names readable."""
    return "".join(
        char if char.isalnum() or char in ("-", "_", ".") else "_"
        for char in str(value)
    )


def resolve_science_path(row: pd.Series, science_root: Path | None) -> Path:
    """Resolve an e2ds file from current metadata and optional Data/science."""
    raw_path = Path(str(row["PATH"]))
    candidates = [raw_path]

    if science_root is not None:
        rel_path = row.get("REL_PATH")
        if pd.notna(rel_path):
            candidates.append(science_root / str(rel_path))
        candidates.append(science_root / str(row.get("NIGHT", "")) / raw_path.name)

    for candidate in candidates:
        if candidate.is_file():
            return candidate

    return raw_path


def read_sources(cube_hdul: fits.HDUList) -> list[str]:
    """Read the Stage 2 SOURCES extension if present."""
    if "SOURCES" not in cube_hdul:
        return []

    source_values = cube_hdul["SOURCES"].data["SOURCE"]
    sources = []
    for value in source_values:
        if isinstance(value, bytes):
            sources.append(value.decode("utf-8").strip())
        else:
            sources.append(str(value).strip())
    return sources


def make_source_hdu(source_names: list[str], original_indices: list[int]) -> fits.BinTableHDU:
    max_len = max([len(name) for name in source_names] + [1])
    source_array = np.array(source_names, dtype=f"S{max_len}")
    return fits.BinTableHDU.from_columns(
        [
            fits.Column(
                name="EXPOSURE",
                format="J",
                array=np.arange(len(source_names), dtype=np.int32),
            ),
            fits.Column(
                name="SOURCE_CUBE_EXPOSURE",
                format="J",
                array=np.asarray(original_indices, dtype=np.int32),
            ),
            fits.Column(name="SOURCE", format=f"{max_len}A", array=source_array),
        ],
        name="SOURCES",
    )


def read_header_metadata(science_path: Path) -> dict[str, object]:
    values = {}
    with fits.open(science_path, memmap=True) as hdul:
        header = hdul[0].header
        for fits_key, column_name in HEADER_KEYS.items():
            values[column_name] = header.get(fits_key)
    return values


def build_parser() -> argparse.ArgumentParser:
    project_dir = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description="Build date-oriented telluric cube products and metadata."
    )
    parser.add_argument(
        "--metadata-final",
        default=str(project_dir / "Output" / "phase1" / "metadata_final.parquet"),
        help="Stage 1 final metadata table.",
    )
    parser.add_argument(
        "--spectra-root",
        default=str(project_dir / "Output" / "spectra"),
        help="Stage 2 object-oriented spectra output root.",
    )
    parser.add_argument(
        "--data-telluric",
        default=str(project_dir / "Data" / "telluric"),
        help="Stage 3 date-oriented telluric dataset root.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(project_dir / "Output" / "stage3"),
        help="Directory for Stage 3 tables and manifest files.",
    )
    parser.add_argument(
        "--science-root",
        default=str(project_dir / "Data" / "science"),
        help="Science root used to resolve REL_PATH values.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing date-oriented cube files.",
    )
    parser.add_argument(
        "--no-csv",
        action="store_true",
        help="Do not write the CSV copy of the Stage 3 index table.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()

    metadata_path = Path(args.metadata_final)
    spectra_root = Path(args.spectra_root)
    data_telluric = Path(args.data_telluric)
    output_dir = Path(args.output_dir)
    science_root = Path(args.science_root) if args.science_root else None

    if not metadata_path.is_file():
        raise FileNotFoundError(f"metadata_final table not found: {metadata_path}")
    if not spectra_root.is_dir():
        raise FileNotFoundError(f"Stage 2 spectra root not found: {spectra_root}")
    if science_root is not None and not science_root.is_dir():
        raise FileNotFoundError(f"Science root not found: {science_root}")

    data_telluric.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_parquet(metadata_path)
    e2ds = metadata[metadata["PRODUCT_TYPE"] == "e2ds_A"].copy()
    e2ds["SOURCE_BASENAME"] = e2ds["PATH"].map(lambda value: Path(str(value)).name)

    index_rows: list[dict[str, object]] = []
    warnings: list[str] = []
    packaged_cubes = 0

    for object_name, object_rows in e2ds.groupby("OBJECT", sort=True):
        object_name = str(object_name)
        safe_object = safe_name(object_name)
        source_rows = {
            row["SOURCE_BASENAME"]: row
            for _, row in object_rows.iterrows()
        }

        source_cube = (
            spectra_root
            / object_name
            / object_name
            / "tell_spec"
            / f"{safe_object}_telluric_cube.fits"
        )
        if not source_cube.is_file():
            warnings.append(f"Missing Stage 2 cube for {object_name}: {source_cube}")
            continue

        with fits.open(source_cube, memmap=True) as hdul:
            cube_data = hdul[0].data
            cube_header = hdul[0].header.copy()
            source_names = read_sources(hdul)

            if not source_names:
                source_names = [
                    row["SOURCE_BASENAME"]
                    for _, row in object_rows.sort_values("SOURCE_BASENAME").iterrows()
                ]

            if len(source_names) != cube_data.shape[0]:
                warnings.append(
                    f"Source count mismatch for {object_name}: "
                    f"{len(source_names)} source rows vs cube axis {cube_data.shape[0]}"
                )
                continue

            exposure_records = []
            for source_cube_index, source_name in enumerate(source_names):
                row = source_rows.get(source_name)
                if row is None:
                    warnings.append(
                        f"No metadata row for {object_name} cube source {source_name}"
                    )
                    continue
                night = str(row["NIGHT"])
                science_path = resolve_science_path(row, science_root)
                exposure_records.append(
                    {
                        "night": night,
                        "source_name": source_name,
                        "source_cube_index": source_cube_index,
                        "metadata_row": row,
                        "science_path": science_path,
                    }
                )

            records_by_night: dict[str, list[dict[str, object]]] = defaultdict(list)
            for record in exposure_records:
                records_by_night[str(record["night"])].append(record)

            for night in sorted(records_by_night):
                night_records = records_by_night[night]
                selected_indices = [
                    int(record["source_cube_index"])
                    for record in night_records
                ]
                packaged_dir = data_telluric / str(night)
                packaged_dir.mkdir(parents=True, exist_ok=True)
                packaged_cube = packaged_dir / f"{safe_object}_{night}_telluric_cube.fits"

                if packaged_cube.exists() and not args.overwrite:
                    raise FileExistsError(
                        f"Stage 3 cube already exists: {packaged_cube}. "
                        "Use --overwrite to replace it."
                    )

                header = cube_header.copy()
                header["OBJECT"] = object_name
                header["NIGHT"] = str(night)
                header["TELLCUBE"] = (True, "Data stored as exposure x order x pixel cube")
                header["NSPECTRA"] = (
                    len(selected_indices),
                    "Number of spectra/exposures in this night cube",
                )
                header["SRC_CUBE"] = (
                    source_cube.name,
                    "Source Stage 2 object cube filename",
                )

                selected_source_names = [
                    str(record["source_name"])
                    for record in night_records
                ]
                fits.HDUList(
                    [
                        fits.PrimaryHDU(
                            np.asarray(cube_data[selected_indices], dtype=np.float32),
                            header=header,
                        ),
                        make_source_hdu(selected_source_names, selected_indices),
                    ]
                ).writeto(packaged_cube, overwrite=True)
                packaged_cubes += 1

                for exposure_index, record in enumerate(night_records):
                    row = record["metadata_row"]
                    science_path = Path(record["science_path"])
                    header_values = read_header_metadata(science_path)
                    rel_cube = packaged_cube.relative_to(data_telluric)
                    index_rows.append(
                        {
                            "OBJECT": object_name,
                            "OBJ_ID": row.get("OBJ_ID"),
                            "NIGHT": str(night),
                            "SOURCE_E2DS": record["source_name"],
                            "SOURCE_REL_PATH": row.get("REL_PATH"),
                            "SOURCE_PATH": str(science_path),
                            "TELLURIC_CUBE_REL_PATH": str(rel_cube),
                            "TELLURIC_CUBE_PATH": str(packaged_cube),
                            "CUBE_EXPOSURE_INDEX": exposure_index,
                            "SOURCE_CUBE_PATH": str(source_cube),
                            "SOURCE_CUBE_EXPOSURE_INDEX": int(
                                record["source_cube_index"]
                            ),
                            "NSPECTRA": len(selected_indices),
                            "NORDERS": int(cube_data.shape[1]),
                            "NPIX": int(cube_data.shape[2]),
                            **header_values,
                        }
                    )

    index = pd.DataFrame(index_rows)
    index_path = output_dir / "telluric_cube_index.parquet"
    index.to_parquet(index_path, index=False)

    csv_path = None
    if not args.no_csv:
        csv_path = output_dir / "telluric_cube_index.csv"
        index.to_csv(csv_path, index=False)

    manifest = {
        "metadata_final": str(metadata_path),
        "spectra_root": str(spectra_root),
        "data_telluric": str(data_telluric),
        "output_dir": str(output_dir),
        "input_e2ds": int(len(e2ds)),
        "packaged_cubes": packaged_cubes,
        "index_rows": int(len(index)),
        "index_path": str(index_path),
        "csv_path": str(csv_path) if csv_path else None,
        "warnings": warnings,
    }
    manifest_path = output_dir / "telluric_stage3_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))

    print(f"Packaged night/object cubes: {packaged_cubes}")
    print(f"Index rows: {len(index)}")
    print(f"Index table: {index_path}")
    if csv_path:
        print(f"CSV table: {csv_path}")
    print(f"Manifest: {manifest_path}")
    if warnings:
        print(f"Warnings: {len(warnings)}")
        for warning in warnings[:10]:
            print(f"WARNING: {warning}")

    return 0 if not warnings else 1


if __name__ == "__main__":
    raise SystemExit(main())
