#!/usr/bin/env python3
"""
Merge benchmark JSON snapshots and produce a consolidated master summary.

Defaults:
- input files: Results/benchmark_*.json (excluding benchmark_latest.json)
- output file: Results/master_benchmark_summary_YYYYMMDD-HHMMSS.json

Merge rule:
- rows are keyed by (Tier, Key, Solver)
- the newest file wins on conflicts
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


Row = Dict[str, object]
RowKey = Tuple[str, str, str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-dir",
        default="Results",
        help="Directory containing benchmark_*.json files (default: Results)",
    )
    parser.add_argument(
        "--inputs",
        nargs="*",
        default=None,
        help=(
            "Explicit benchmark JSON files to merge. "
            "If omitted, auto-discovers benchmark_*.json in --results-dir."
        ),
    )
    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Output JSON path. "
            "Default: Results/master_benchmark_summary_YYYYMMDD-HHMMSS.json"
        ),
    )
    parser.add_argument(
        "--print-files",
        action="store_true",
        help="Print the input files selected for merge.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail on unreadable/invalid JSON input files instead of skipping.",
    )
    return parser.parse_args()


def discover_input_files(results_dir: Path) -> List[Path]:
    files = sorted(results_dir.glob("benchmark_*.json"))
    return [f for f in files if f.name != "benchmark_latest.json"]


def load_rows(path: Path) -> List[Row]:
    data = json.loads(path.read_text())
    if not isinstance(data, list):
        raise ValueError(f"{path} is not a benchmark row list JSON.")
    rows: List[Row] = []
    for i, row in enumerate(data):
        if not isinstance(row, dict):
            raise ValueError(f"{path} row #{i} is not an object.")
        rows.append(row)
    return rows


def row_key(row: Row) -> RowKey:
    tier = str(row.get("Tier", ""))
    key = str(row.get("Key", ""))
    solver = str(row.get("Solver", ""))
    return tier, key, solver


def merge_rows(files: Iterable[Path], strict: bool = False) -> Tuple[List[Row], List[str], List[Dict[str, str]]]:
    # Newest file should win. We process oldest -> newest and overwrite by key.
    merged: Dict[RowKey, Row] = {}
    file_names: List[str] = []
    skipped: List[Dict[str, str]] = []
    for path in files:
        try:
            if not path.exists():
                raise FileNotFoundError(str(path))
            if path.stat().st_size == 0:
                raise ValueError("empty file")
            rows = load_rows(path)
        except Exception as exc:  # pylint: disable=broad-exception-caught
            if strict:
                raise
            skipped.append({"file": path.name, "reason": str(exc)})
            continue

        file_names.append(path.name)
        for row in rows:
            merged[row_key(row)] = row
    rows = list(merged.values())
    rows.sort(key=lambda r: (str(r.get("Tier", "")), str(r.get("Key", "")), str(r.get("Solver", ""))))
    return rows, file_names, skipped


def summarize(rows: List[Row], files_used: List[str], skipped_files: List[Dict[str, str]]) -> Dict[str, object]:
    status_counts = Counter(str(r.get("Status", "UNKNOWN")) for r in rows)

    by_solver: Dict[str, Counter] = defaultdict(Counter)
    for r in rows:
        by_solver[str(r.get("Solver", "UNKNOWN"))][str(r.get("Status", "UNKNOWN"))] += 1

    by_tier_rows: Dict[str, List[Row]] = defaultdict(list)
    for r in rows:
        by_tier_rows[str(r.get("Tier", "unknown"))].append(r)

    by_tier = {}
    for tier, tier_rows in sorted(by_tier_rows.items()):
        case_count = len({str(r.get("Key", "")) for r in tier_rows})
        tier_status = Counter(str(r.get("Status", "UNKNOWN")) for r in tier_rows)
        tier_solver: Dict[str, Counter] = defaultdict(Counter)
        for r in tier_rows:
            tier_solver[str(r.get("Solver", "UNKNOWN"))][str(r.get("Status", "UNKNOWN"))] += 1

        by_tier[tier] = {
            "cases": case_count,
            "rows": len(tier_rows),
            "status": dict(tier_status),
            "by_solver": {solver: dict(counts) for solver, counts in sorted(tier_solver.items())},
        }

    non_ok = [r for r in rows if str(r.get("Status", "UNKNOWN")) != "OK"]

    summary = {
        "files": files_used,
        "skipped_files": skipped_files,
        "master_cases": len({str(r.get("Key", "")) for r in rows}),
        "master_rows": len(rows),
        "status": dict(status_counts),
        "by_solver": {solver: dict(counts) for solver, counts in sorted(by_solver.items())},
        "by_tier": by_tier,
        "non_ok": non_ok,
    }
    return summary


def main() -> int:
    args = parse_args()
    results_dir = Path(args.results_dir).resolve()
    if not results_dir.exists():
        raise SystemExit(f"Results dir not found: {results_dir}")

    if args.inputs:
        files = [Path(p).resolve() for p in args.inputs]
    else:
        files = discover_input_files(results_dir)

    if not files:
        raise SystemExit("No benchmark JSON files found to merge.")

    if args.print_files:
        print("Merging files:")
        for f in files:
            print(f"  - {f}")

    rows, files_used, skipped_files = merge_rows(files, strict=args.strict)
    summary = summarize(rows, files_used, skipped_files)

    if args.output:
        out_path = Path(args.output).resolve()
    else:
        ts = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
        out_path = results_dir / f"master_benchmark_summary_{ts}.json"

    out_path.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n")
    print(f"Wrote summary: {out_path}")
    print(f"Rows: {summary['master_rows']} | Cases: {summary['master_cases']}")
    print(f"Status counts: {summary['status']}")
    if skipped_files:
        print(f"Skipped files: {len(skipped_files)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
