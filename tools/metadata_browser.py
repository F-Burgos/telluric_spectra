#!/usr/bin/env python3
"""Small browser-based metadata table explorer."""

from __future__ import annotations

import argparse
import html
import os
import re
import socketserver
import sys
import webbrowser
from difflib import SequenceMatcher
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler
from pathlib import Path
from urllib.parse import parse_qs, urlencode, urlparse

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_METADATA = PROJECT_DIR / "Output" / "phase1" / "metadata_final.parquet"
DEFAULT_HOST = "127.0.0.1"
DEFAULT_PORT = 8765
DEFAULT_LIMIT = 200


def normalize_token(value: object) -> str:
    """Normalize object labels for punctuation-insensitive matching."""
    return "".join(char for char in str(value).casefold() if char.isalnum())


def normalized_match(label: object, query_token: str) -> bool:
    """Return True for punctuation-insensitive and minor-spelling object matches."""
    label_token = normalize_token(label)
    if not query_token:
        return True
    if query_token in label_token or label_token in query_token:
        return True
    if len(query_token) < 5 or len(label_token) < 5:
        return False
    return SequenceMatcher(None, query_token, label_token).ratio() >= 0.88


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Serve an interactive browser UI for filtering metadata parquet tables."
    )
    parser.add_argument(
        "--metadata",
        default=str(DEFAULT_METADATA),
        help="Metadata parquet table. Default: Output/phase1/metadata_final.parquet",
    )
    parser.add_argument(
        "--host",
        default=DEFAULT_HOST,
        help=f"Host/interface to bind. Default: {DEFAULT_HOST}",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=DEFAULT_PORT,
        help=f"Port to bind. Default: {DEFAULT_PORT}",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Do not automatically open the browser.",
    )
    return parser.parse_args()


class MetadataStore:
    """Lazy metadata loader with simple mtime-based caching."""

    def __init__(self, metadata_path: Path):
        self.metadata_path = metadata_path
        self._mtime: float | None = None
        self._dataframe: pd.DataFrame | None = None

    def load(self) -> pd.DataFrame:
        if not self.metadata_path.is_file():
            raise FileNotFoundError(f"Metadata table not found: {self.metadata_path}")

        mtime = self.metadata_path.stat().st_mtime
        if self._dataframe is None or self._mtime != mtime:
            self._dataframe = pd.read_parquet(self.metadata_path)
            self._mtime = mtime
        return self._dataframe.copy()


def first_value(params: dict[str, list[str]], name: str, default: str = "") -> str:
    return params.get(name, [default])[0].strip()


def sorted_unique(df: pd.DataFrame, column: str) -> list[str]:
    if column not in df.columns:
        return []
    values = df[column].dropna().astype(str).drop_duplicates().sort_values()
    return values.tolist()


def filter_metadata(df: pd.DataFrame, params: dict[str, list[str]]) -> tuple[pd.DataFrame, list[str]]:
    """Apply UI filters and return filtered rows plus human-readable notes."""
    filtered = df.copy()
    notes = []

    object_query = first_value(params, "star")
    object_column = first_value(params, "object_column", "OBJECT") or "OBJECT"
    match_mode = first_value(params, "match_mode", "normalized")
    product_type = first_value(params, "product_type")
    night = first_value(params, "night")
    date_from = first_value(params, "date_from")
    date_to = first_value(params, "date_to")

    if object_query and object_column in filtered.columns:
        labels = filtered[object_column].astype(str)
        if match_mode == "exact":
            mask = labels.str.casefold().eq(object_query.casefold())
        elif match_mode == "contains":
            mask = labels.str.contains(object_query, case=False, regex=False, na=False)
        else:
            query_token = normalize_token(object_query)
            mask = labels.map(lambda label: normalized_match(label, query_token))
        filtered = filtered.loc[mask]
        notes.append(f"{object_column} matches '{object_query}' ({match_mode})")

    if product_type and "PRODUCT_TYPE" in filtered.columns:
        filtered = filtered.loc[filtered["PRODUCT_TYPE"].astype(str) == product_type]
        notes.append(f"PRODUCT_TYPE = {product_type}")

    if night and "NIGHT" in filtered.columns:
        filtered = filtered.loc[filtered["NIGHT"].astype(str) == night]
        notes.append(f"NIGHT = {night}")

    if (date_from or date_to) and "NIGHT" in filtered.columns:
        dates = pd.to_datetime(filtered["NIGHT"].astype(str), errors="coerce")
        mask = dates.notna()
        if date_from:
            start = pd.to_datetime(date_from)
            mask &= dates >= start
            notes.append(f"NIGHT >= {date_from}")
        if date_to:
            end = pd.to_datetime(date_to) + pd.Timedelta(days=1)
            mask &= dates < end
            notes.append(f"NIGHT <= {date_to}")
        filtered = filtered.loc[mask]

    return filtered, notes


def numeric_param(params: dict[str, list[str]], name: str, default: int) -> int:
    raw = first_value(params, name, str(default))
    try:
        return max(1, int(raw))
    except ValueError:
        return default


def render_select(name: str, values: list[str], selected: str, placeholder: str) -> str:
    options = [f'<option value="">{html.escape(placeholder)}</option>']
    for value in values:
        escaped = html.escape(value)
        is_selected = " selected" if value == selected else ""
        options.append(f'<option value="{escaped}"{is_selected}>{escaped}</option>')
    return f'<select name="{html.escape(name)}">{"".join(options)}</select>'


def format_stats(df: pd.DataFrame) -> list[str]:
    stats = [f"Rows: {len(df)}"]
    for column in ["OBJECT", "OBJ_ID", "NIGHT", "PRODUCT_TYPE"]:
        if column in df.columns:
            stats.append(f"Unique {column}: {df[column].nunique(dropna=True)}")
    if "NIGHT" in df.columns:
        nights = pd.to_datetime(df["NIGHT"].astype(str), errors="coerce").dropna()
        if not nights.empty:
            stats.append(f"Night span: {nights.min().date()} to {nights.max().date()}")
    return stats


def dataframe_to_html(df: pd.DataFrame) -> str:
    if df.empty:
        return "<p>No rows matched the current query.</p>"
    return df.to_html(
        index=False,
        escape=True,
        classes="metadata-table",
        border=0,
        max_cols=None,
    )


def render_page(
    metadata_path: Path,
    df: pd.DataFrame,
    filtered: pd.DataFrame,
    shown: pd.DataFrame,
    params: dict[str, list[str]],
    notes: list[str],
    page: int,
    limit: int,
) -> bytes:
    star = html.escape(first_value(params, "star"))
    object_column = first_value(params, "object_column", "OBJECT") or "OBJECT"
    match_mode = first_value(params, "match_mode", "normalized") or "normalized"
    product_type = first_value(params, "product_type")
    night = first_value(params, "night")
    date_from = html.escape(first_value(params, "date_from"))
    date_to = html.escape(first_value(params, "date_to"))

    columns = [column for column in ["OBJECT", "OBJ_ID"] if column in df.columns]
    if not columns:
        columns = df.columns[:1].tolist()

    object_column_options = []
    for column in columns:
        selected = " selected" if column == object_column else ""
        object_column_options.append(
            f'<option value="{html.escape(column)}"{selected}>{html.escape(column)}</option>'
        )

    mode_options = []
    for value, label in [
        ("normalized", "Normalized contains"),
        ("contains", "Text contains"),
        ("exact", "Exact"),
    ]:
        selected = " selected" if value == match_mode else ""
        mode_options.append(f'<option value="{value}"{selected}>{label}</option>')

    product_select = render_select(
        "product_type",
        sorted_unique(df, "PRODUCT_TYPE"),
        product_type,
        "Any product",
    )
    night_select = render_select("night", sorted_unique(df, "NIGHT"), night, "Any night")

    csv_params = {key: value[0] for key, value in params.items() if value and value[0]}
    csv_params["format"] = "csv"
    csv_link = f"/?{urlencode(csv_params)}"

    previous_page = max(1, page - 1)
    next_page = page + 1
    base_params = {key: value[0] for key, value in params.items() if value and value[0]}
    base_params["limit"] = str(limit)
    previous_params = dict(base_params, page=str(previous_page))
    next_params = dict(base_params, page=str(next_page))

    stats_html = "".join(f"<li>{html.escape(stat)}</li>" for stat in format_stats(filtered))
    notes_html = "".join(f"<li>{html.escape(note)}</li>" for note in notes)

    body = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Telluric metadata browser</title>
  <style>
    body {{
      font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      margin: 24px;
      color: #17202a;
      background: #fbfcfc;
    }}
    form {{
      display: grid;
      grid-template-columns: repeat(4, minmax(180px, 1fr));
      gap: 12px;
      align-items: end;
      padding: 16px;
      background: #eef5ff;
      border: 1px solid #c9d8f2;
      border-radius: 8px;
      margin-bottom: 16px;
    }}
    label {{
      display: flex;
      flex-direction: column;
      gap: 4px;
      font-size: 0.9rem;
      font-weight: 600;
    }}
    input, select, button {{
      font: inherit;
      padding: 7px 8px;
      border-radius: 6px;
      border: 1px solid #aeb6bf;
      background: white;
    }}
    button {{
      background: #1f77b4;
      color: white;
      border-color: #1f77b4;
      cursor: pointer;
      font-weight: 700;
    }}
    .toolbar {{
      display: flex;
      gap: 12px;
      flex-wrap: wrap;
      align-items: center;
      margin: 12px 0;
    }}
    .panel {{
      padding: 12px 16px;
      border: 1px solid #d5dbdb;
      border-radius: 8px;
      background: white;
      margin-bottom: 16px;
    }}
    .metadata-table {{
      border-collapse: collapse;
      width: 100%;
      font-size: 0.85rem;
      background: white;
    }}
    .metadata-table th, .metadata-table td {{
      border: 1px solid #d5dbdb;
      padding: 5px 7px;
      vertical-align: top;
      white-space: nowrap;
    }}
    .metadata-table th {{
      position: sticky;
      top: 0;
      background: #ecf0f1;
      z-index: 1;
    }}
    .table-wrap {{
      overflow: auto;
      max-height: 70vh;
      border: 1px solid #d5dbdb;
      border-radius: 8px;
    }}
    code {{ background: #f4f6f7; padding: 2px 4px; border-radius: 4px; }}
    a {{ color: #1f618d; }}
  </style>
</head>
<body>
  <h1>Telluric metadata browser</h1>
  <p>Table: <code>{html.escape(str(metadata_path))}</code></p>

  <form method="get" action="/">
    <label>
      Star/Object query
      <input name="star" value="{star}" placeholder="HD10700, object alias, ..." autofocus>
    </label>
    <label>
      Object column
      <select name="object_column">{"".join(object_column_options)}</select>
    </label>
    <label>
      Match mode
      <select name="match_mode">{"".join(mode_options)}</select>
    </label>
    <label>
      Product type
      {product_select}
    </label>
    <label>
      Exact night
      {night_select}
    </label>
    <label>
      Night/date from
      <input name="date_from" value="{date_from}" placeholder="2015-10-02">
    </label>
    <label>
      Night/date to
      <input name="date_to" value="{date_to}" placeholder="2015-10-04">
    </label>
    <label>
      Rows per page
      <input name="limit" value="{limit}" type="number" min="1">
    </label>
    <button type="submit">Query metadata</button>
  </form>

  <div class="panel">
    <h2>Result summary</h2>
    <ul>{stats_html}</ul>
    <p>Showing rows {(page - 1) * limit + 1 if len(filtered) else 0}
       to {min(page * limit, len(filtered))} of {len(filtered)}.</p>
    <ul>{notes_html or "<li>No filters active.</li>"}</ul>
  </div>

  <div class="toolbar">
    <a href="/?{urlencode(previous_params)}">Previous page</a>
    <span>Page {page}</span>
    <a href="/?{urlencode(next_params)}">Next page</a>
    <a href="{csv_link}">Download filtered CSV</a>
    <a href="/">Reset</a>
  </div>

  <div class="table-wrap">
    {dataframe_to_html(shown)}
  </div>
</body>
</html>
"""
    return body.encode("utf-8")


class MetadataRequestHandler(BaseHTTPRequestHandler):
    store: MetadataStore

    def log_message(self, format: str, *args: object) -> None:
        sys.stderr.write(f"[metadata-browser] {format % args}\n")

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        if parsed.path not in {"/", "/index.html"}:
            self.send_error(HTTPStatus.NOT_FOUND)
            return

        params = parse_qs(parsed.query)
        try:
            df = self.store.load()
            filtered, notes = filter_metadata(df, params)
            page = numeric_param(params, "page", 1)
            limit = numeric_param(params, "limit", DEFAULT_LIMIT)
            start = (page - 1) * limit
            shown = filtered.iloc[start : start + limit]

            if first_value(params, "format") == "csv":
                payload = filtered.to_csv(index=False).encode("utf-8")
                self.send_response(HTTPStatus.OK)
                self.send_header("Content-Type", "text/csv; charset=utf-8")
                self.send_header(
                    "Content-Disposition",
                    'attachment; filename="metadata_query.csv"',
                )
                self.send_header("Content-Length", str(len(payload)))
                self.end_headers()
                self.wfile.write(payload)
                return

            payload = render_page(
                metadata_path=self.store.metadata_path,
                df=df,
                filtered=filtered,
                shown=shown,
                params=params,
                notes=notes,
                page=page,
                limit=limit,
            )
            self.send_response(HTTPStatus.OK)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)
        except Exception as exc:
            payload = (
                "<!doctype html><html><body>"
                "<h1>Metadata browser error</h1>"
                f"<pre>{html.escape(str(exc))}</pre>"
                "</body></html>"
            ).encode("utf-8")
            self.send_response(HTTPStatus.INTERNAL_SERVER_ERROR)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)


def main() -> int:
    args = parse_args()
    metadata_path = Path(args.metadata)
    MetadataRequestHandler.store = MetadataStore(metadata_path)

    class ReusableTCPServer(socketserver.TCPServer):
        allow_reuse_address = True

    with ReusableTCPServer((args.host, args.port), MetadataRequestHandler) as server:
        url = f"http://{args.host}:{args.port}/"
        print(f"Serving metadata browser at {url}")
        print(f"Metadata table: {metadata_path}")
        print("Press Ctrl+C to stop.")
        if not args.no_browser:
            webbrowser.open(url)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            print("\nStopping metadata browser.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
