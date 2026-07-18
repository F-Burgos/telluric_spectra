#!/usr/bin/env python3
"""Plot sky positions for one or more OBJECT groups from a metadata table."""

from __future__ import annotations

import argparse
import os
import re
import textwrap
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_MPLCONFIGDIR = PROJECT_DIR / "Output" / ".matplotlib"
DEFAULT_MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(DEFAULT_MPLCONFIGDIR))

import matplotlib.pyplot as plt


AUTO_DATE_COLUMNS = ("NIGHT", "DATE_OBS", "DATE-OBS", "MJD_OBS", "MJD-OBS")


def parse_objects(values: list[str]) -> list[str]:
    """Accept repeated values and comma-separated object lists."""
    objects: list[str] = []
    for value in values:
        objects.extend(part.strip() for part in value.split(",") if part.strip())
    return objects


def safe_name(value: str) -> str:
    """Return a filesystem-safe label for plot filenames."""
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return cleaned.strip("_") or "object"


def match_objects(
    metadata: pd.DataFrame,
    requested_objects: list[str],
    object_column: str,
    contains: bool,
) -> pd.DataFrame:
    """Filter the metadata table by OBJECT-like labels."""
    labels = metadata[object_column].astype(str)
    mask = pd.Series(False, index=metadata.index)

    for object_name in requested_objects:
        if contains:
            mask |= labels.str.contains(object_name, case=False, regex=False, na=False)
        else:
            mask |= labels.str.casefold().eq(object_name.casefold())

    return metadata.loc[mask].copy()


def available_object_hint(metadata: pd.DataFrame, object_column: str, query: str) -> str:
    """Return a short hint with available object names similar to the query."""
    labels = metadata[object_column].dropna().astype(str).drop_duplicates().sort_values()
    query_folded = query.casefold()
    candidates = [
        label
        for label in labels
        if query_folded in label.casefold() or label.casefold() in query_folded
    ]
    if not candidates:
        candidates = labels.head(10).tolist()
    return ", ".join(candidates[:10])


def resolve_date_column(metadata: pd.DataFrame, requested_column: str) -> str | None:
    """Choose a date column for filtering/statistics."""
    if requested_column != "auto":
        if requested_column not in metadata.columns:
            raise ValueError(f"Requested date column not found: {requested_column}")
        return requested_column

    for column in AUTO_DATE_COLUMNS:
        if column in metadata.columns:
            return column
    return None


def add_observation_dates(metadata: pd.DataFrame, date_column: str | None) -> pd.DataFrame:
    """Attach a normalized helper datetime column when a date-like column exists."""
    metadata = metadata.copy()
    if date_column is None:
        metadata["_OBS_DATE"] = pd.NaT
        return metadata

    if date_column in {"MJD_OBS", "MJD-OBS"}:
        mjd = pd.to_numeric(metadata[date_column], errors="coerce")
        metadata["_OBS_DATE"] = pd.to_datetime(
            mjd,
            unit="D",
            origin="1858-11-17",
            errors="coerce",
        )
    else:
        metadata["_OBS_DATE"] = pd.to_datetime(
            metadata[date_column],
            errors="coerce",
            utc=False,
        )
    return metadata


def filter_by_date(
    metadata: pd.DataFrame,
    date_from: str | None,
    date_to: str | None,
    date_column: str | None,
) -> pd.DataFrame:
    """Filter rows by inclusive observation date limits."""
    if not date_from and not date_to:
        return metadata

    if date_column is None:
        raise ValueError(
            "Date filtering was requested, but no date column was found. "
            "Use --date-column to select one explicitly."
        )

    filtered = metadata.copy()
    if date_from:
        start = pd.to_datetime(date_from)
        filtered = filtered.loc[filtered["_OBS_DATE"] >= start]
    if date_to:
        end = pd.to_datetime(date_to) + pd.Timedelta(days=1)
        filtered = filtered.loc[filtered["_OBS_DATE"] < end]
    return filtered


def format_date_span(rows: pd.DataFrame) -> str:
    """Return a compact observation date span for stats text."""
    dates = rows["_OBS_DATE"].dropna() if "_OBS_DATE" in rows.columns else pd.Series([])
    if dates.empty:
        return "date span unavailable"
    return f"{dates.min().date()} to {dates.max().date()}"


def pluralize(count: int, singular: str, plural: str | None = None) -> str:
    """Return a count with a simple English singular/plural label."""
    if count == 1:
        return f"{count} {singular}"
    return f"{count} {plural or singular + 's'}"


def add_side_statistics(
    info_ax,
    group: pd.DataFrame,
    object_column: str,
    date_column: str | None,
) -> None:
    """Add observation statistics in a dedicated side panel."""
    lines = [
        "Statistics",
        "",
        f"Total observations: {len(group)}",
        f"Date column: {date_column or 'unavailable'}",
        f"Date span: {format_date_span(group)}",
    ]

    if "NIGHT" in group.columns:
        lines.append(f"Unique nights: {group['NIGHT'].nunique()}")
    if "OBJ_ID" in group.columns:
        lines.append(f"Unique OBJ_ID: {group['OBJ_ID'].nunique()}")

    lines.extend(["", "Objects"])
    for object_name, object_rows in group.groupby(object_column, sort=True):
        object_lines = [f"{object_name}:"]
        object_lines.append(f"  observations: {len(object_rows)}")
        if "NIGHT" in object_rows.columns:
            object_lines.append(f"  nights: {object_rows['NIGHT'].nunique()}")
        if "OBJ_ID" in object_rows.columns:
            object_lines.append(f"  OBJ_ID values: {object_rows['OBJ_ID'].nunique()}")
        object_lines.append(f"  date span: {format_date_span(object_rows)}")
        object_lines.append(
            f"  RA: {object_rows['RA'].min():.6f} to {object_rows['RA'].max():.6f}"
        )
        object_lines.append(
            f"  DEC: {object_rows['DEC'].min():.6f} to {object_rows['DEC'].max():.6f}"
        )
        lines.extend(object_lines)

    wrapped_lines = []
    for line in lines:
        if line.startswith("  "):
            wrapped_lines.extend(
                textwrap.wrap(line, width=34, subsequent_indent="    ") or [line]
            )
        else:
            wrapped_lines.extend(textwrap.wrap(line, width=34) or [line])

    info_ax.axis("off")
    info_ax.text(
        0.0,
        0.98,
        "\n".join(wrapped_lines),
        transform=info_ax.transAxes,
        va="top",
        ha="left",
        fontsize=9.5,
        linespacing=1.25,
    )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Plot RA/DEC sky positions for one or more canonical OBJECT groups "
            "from a metadata parquet table."
        )
    )
    parser.add_argument(
        "--metadata",
        default=str(PROJECT_DIR / "Output" / "phase1" / "metadata_final.parquet"),
        help="Metadata parquet table. Default: Output/phase1/metadata_final.parquet",
    )
    parser.add_argument(
        "--object",
        action="append",
        required=True,
        help=(
            "Canonical OBJECT to plot. Can be repeated or comma-separated, "
            'for example: --object Barnards_Star --object HD20794 or --object "HD20794,HD39091".'
        ),
    )
    parser.add_argument(
        "--object-column",
        default="OBJECT",
        help="Column containing the object label to filter. Default: OBJECT",
    )
    parser.add_argument(
        "--contains",
        action="store_true",
        help="Use case-insensitive substring matching instead of exact matching.",
    )
    parser.add_argument(
        "--date-column",
        default="auto",
        help=(
            "Date column used for filtering/statistics. Default: auto "
            f"({', '.join(AUTO_DATE_COLUMNS)})."
        ),
    )
    parser.add_argument(
        "--date-from",
        default=None,
        help="Inclusive first observation date to plot, for example 2012-01-01.",
    )
    parser.add_argument(
        "--date-to",
        default=None,
        help="Inclusive last observation date to plot, for example 2012-12-31.",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output PNG path. Default: Output/figures/sky_positions_<objects>.png",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Optional custom plot title.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Open an interactive matplotlib window instead of only saving the PNG.",
    )
    parser.add_argument(
        "--no-invert-ra",
        action="store_true",
        help="Do not invert the RA axis. By default RA increases toward the left.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.7,
        help="Scatter marker transparency. Default: 0.7",
    )
    parser.add_argument(
        "--marker-size",
        type=float,
        default=24.0,
        help="Scatter marker size. Default: 24",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    metadata_path = Path(args.metadata)
    requested_objects = parse_objects(args.object)

    if not metadata_path.is_file():
        raise FileNotFoundError(f"Metadata table not found: {metadata_path}")

    metadata = pd.read_parquet(metadata_path)
    required_columns = {"RA", "DEC", args.object_column}
    missing_columns = sorted(required_columns.difference(metadata.columns))
    if missing_columns:
        raise ValueError(
            f"Metadata table is missing required columns: {', '.join(missing_columns)}"
        )
    date_column = resolve_date_column(metadata, args.date_column)
    metadata = add_observation_dates(metadata, date_column)

    selected = match_objects(
        metadata=metadata,
        requested_objects=requested_objects,
        object_column=args.object_column,
        contains=args.contains,
    )

    if selected.empty:
        hints = [
            f"{name}: {available_object_hint(metadata, args.object_column, name)}"
            for name in requested_objects
        ]
        raise ValueError(
            "No rows matched the requested object(s). Similar/available labels: "
            + " | ".join(hints)
        )

    selected = filter_by_date(
        metadata=selected,
        date_from=args.date_from,
        date_to=args.date_to,
        date_column=date_column,
    )
    if selected.empty:
        raise ValueError(
            "The requested object(s) were found, but no rows remain after the "
            f"date filter ({args.date_from or 'beginning'} to {args.date_to or 'end'})."
        )

    selected["RA"] = pd.to_numeric(selected["RA"], errors="coerce")
    selected["DEC"] = pd.to_numeric(selected["DEC"], errors="coerce")
    selected = selected.dropna(subset=["RA", "DEC"])
    if selected.empty:
        raise ValueError("Matched rows do not contain finite RA/DEC values.")

    if args.output:
        output_path = Path(args.output)
    else:
        label = "_".join(safe_name(name) for name in requested_objects)
        output_path = PROJECT_DIR / "Output" / "figures" / f"sky_positions_{label}.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, (ax, info_ax) = plt.subplots(
        1,
        2,
        figsize=(12, 6),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [3.0, 1.35]},
    )

    for object_name, object_rows in selected.groupby(args.object_column, sort=True):
        ax.scatter(
            object_rows["RA"],
            object_rows["DEC"],
            s=args.marker_size,
            alpha=args.alpha,
            label=f"{object_name} ({len(object_rows)})",
            edgecolors="none",
        )

    ax.set_xlabel("RA [deg]")
    ax.set_ylabel("DEC [deg]")
    title = args.title or f"Sky positions from {metadata_path.name}: {', '.join(requested_objects)}"
    if args.date_from or args.date_to:
        title += f" ({args.date_from or 'beginning'} to {args.date_to or 'end'})"
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=9)
    add_side_statistics(info_ax, selected, args.object_column, date_column)

    if not args.no_invert_ra:
        ax.invert_xaxis()

    fig.savefig(output_path, dpi=180)
    print(f"Wrote sky-position plot: {output_path}")
    print(f"Matched observations: {len(selected)}")
    print(f"Date column: {date_column or 'unavailable'}")
    print(f"Date span: {format_date_span(selected)}")
    if "NIGHT" in selected.columns:
        print(f"Unique nights: {selected['NIGHT'].nunique()}")
    if "OBJ_ID" in selected.columns:
        print(f"Unique OBJ_ID values: {selected['OBJ_ID'].nunique()}")
    if args.date_from or args.date_to:
        print(
            "Date filter: "
            f"{args.date_from or 'beginning'} to {args.date_to or 'end'}"
        )
    print(
        "Matched OBJECT groups: "
        + ", ".join(
            f"{name} ({len(rows)})"
            for name, rows in selected.groupby(args.object_column, sort=True)
        )
    )
    print("Per-object statistics:")
    for object_name, object_rows in selected.groupby(args.object_column, sort=True):
        line = f"  - {object_name}: {pluralize(len(object_rows), 'observation')}"
        if "NIGHT" in object_rows.columns:
            line += f", {pluralize(object_rows['NIGHT'].nunique(), 'night')}"
        if "OBJ_ID" in object_rows.columns:
            line += (
                f", {pluralize(object_rows['OBJ_ID'].nunique(), 'original OBJ_ID value')}"
            )
        line += f", date span {format_date_span(object_rows)}"
        line += (
            f", RA {object_rows['RA'].min():.6f} to {object_rows['RA'].max():.6f}"
            f", DEC {object_rows['DEC'].min():.6f} to {object_rows['DEC'].max():.6f}"
        )
        print(line)

    if args.show:
        plt.show()
    else:
        plt.close(fig)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
