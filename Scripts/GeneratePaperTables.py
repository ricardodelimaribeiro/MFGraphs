#!/usr/bin/env python3
"""
Generate publication-ready side-by-side benchmark tables from summary JSON files.

Expected input schema matches Results/master_benchmark_summary_*.json:
  - by_tier[tier].by_solver[solver][status] counts
  - status values like OK, TIMEOUT, FAILED

Outputs:
  - Markdown table
  - LaTeX table
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional


DEFAULT_TIERS = ["small", "core", "large", "vlarge", "stress", "paper"]
DEFAULT_SOLVERS = ["CriticalCongestion", "NonLinearSolver", "Monotone"]
STATUS_ORDER = ["OK", "TIMEOUT", "FAILED", "SKIPPED", "UNKNOWN"]

TIER_LABELS = {
    "small": "Small",
    "core": "Core",
    "large": "Large",
    "vlarge": "Very Large",
    "stress": "Stress",
    "paper": "Paper",
    "inconsistent-switching": "Inconsistent Switching",
}

SOLVER_LABELS = {
    "CriticalCongestion": "Critical",
    "NonLinearSolver": "NonLinear",
    "Monotone": "Monotone",
    "DataToEquations": "DataToEquations",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--challenge",
        required=True,
        help="Path to challenge summary JSON (e.g., strict 120s).",
    )
    parser.add_argument(
        "--extended",
        required=True,
        help="Path to extended-budget summary JSON.",
    )
    parser.add_argument(
        "--format",
        choices=["markdown", "latex", "both"],
        default="markdown",
        help="Output format (default: markdown).",
    )
    parser.add_argument(
        "--tiers",
        default=",".join(DEFAULT_TIERS),
        help=f"Comma-separated tiers (default: {','.join(DEFAULT_TIERS)}).",
    )
    parser.add_argument(
        "--solvers",
        default=",".join(DEFAULT_SOLVERS),
        help=f"Comma-separated solvers (default: {','.join(DEFAULT_SOLVERS)}).",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Optional output file path. If omitted, prints to stdout.",
    )
    parser.add_argument(
        "--show-notes",
        action="store_true",
        help="Append notes (e.g., inconsistent-switching compiler validation summary).",
    )
    return parser.parse_args()


def load_summary(path: str) -> Dict[str, object]:
    p = Path(path).resolve()
    data = json.loads(p.read_text())
    if not isinstance(data, dict):
        raise ValueError(f"{p} is not a summary object.")
    return data


def normalize_csv_list(raw: str) -> List[str]:
    return [item.strip() for item in raw.split(",") if item.strip()]


def status_counts(summary: Dict[str, object], tier: str, solver: str) -> Optional[Dict[str, int]]:
    by_tier = summary.get("by_tier", {})
    if not isinstance(by_tier, dict):
        return None
    tier_obj = by_tier.get(tier)
    if not isinstance(tier_obj, dict):
        return None
    by_solver = tier_obj.get("by_solver", {})
    if not isinstance(by_solver, dict):
        return None
    solver_obj = by_solver.get(solver)
    if not isinstance(solver_obj, dict):
        return None
    out: Dict[str, int] = {}
    for k, v in solver_obj.items():
        try:
            out[str(k)] = int(v)
        except Exception:  # pylint: disable=broad-exception-caught
            continue
    return out


def format_counts(counts: Optional[Dict[str, int]]) -> str:
    if counts is None:
        return "—"
    parts: List[str] = []
    seen = set()
    for status in STATUS_ORDER:
        if counts.get(status, 0) > 0:
            parts.append(f"{counts[status]} {status}")
            seen.add(status)
    for status in sorted(counts.keys()):
        if status not in seen and counts.get(status, 0) > 0:
            parts.append(f"{counts[status]} {status}")
    return ", ".join(parts) if parts else "0"


def tier_label(tier: str) -> str:
    return TIER_LABELS.get(tier, tier)


def solver_label(solver: str) -> str:
    return SOLVER_LABELS.get(solver, solver)


def build_markdown(challenge: Dict[str, object], extended: Dict[str, object], tiers: List[str], solvers: List[str]) -> str:
    lines = [
        "| Tier | Solver | 120s Challenge | Extended Budget |",
        "| :--- | :--- | :--- | :--- |",
    ]
    for tier in tiers:
        for solver in solvers:
            c = format_counts(status_counts(challenge, tier, solver))
            e = format_counts(status_counts(extended, tier, solver))
            lines.append(f"| {tier_label(tier)} | {solver_label(solver)} | {c} | {e} |")
    return "\n".join(lines)


def latex_escape(text: str) -> str:
    repl = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    out = text
    for k, v in repl.items():
        out = out.replace(k, v)
    return out


def build_latex(challenge: Dict[str, object], extended: Dict[str, object], tiers: List[str], solvers: List[str]) -> str:
    lines = [
        r"\begin{tabular}{llll}",
        r"\hline",
        r"Tier & Solver & 120s Challenge & Extended Budget \\",
        r"\hline",
    ]
    for tier in tiers:
        for solver in solvers:
            c = latex_escape(format_counts(status_counts(challenge, tier, solver)))
            e = latex_escape(format_counts(status_counts(extended, tier, solver)))
            lines.append(
                f"{latex_escape(tier_label(tier))} & {latex_escape(solver_label(solver))} & {c} & {e} \\\\"
            )
    lines.extend([r"\hline", r"\end{tabular}"])
    return "\n".join(lines)


def build_notes(challenge: Dict[str, object], extended: Dict[str, object]) -> str:
    tier = "inconsistent-switching"
    solver = "DataToEquations"
    c = format_counts(status_counts(challenge, tier, solver))
    e = format_counts(status_counts(extended, tier, solver))
    return (
        "Notes:\n"
        f"- Compiler validation tier ({tier_label(tier)}) via {solver_label(solver)}: "
        f"challenge={c}; extended={e}.\n"
        "- Interpretation: these are expected preprocessing rejections by design."
    )


def main() -> int:
    args = parse_args()
    challenge = load_summary(args.challenge)
    extended = load_summary(args.extended)
    tiers = normalize_csv_list(args.tiers)
    solvers = normalize_csv_list(args.solvers)

    chunks: List[str] = []
    if args.format in {"markdown", "both"}:
        chunks.append(build_markdown(challenge, extended, tiers, solvers))
    if args.format in {"latex", "both"}:
        chunks.append(build_latex(challenge, extended, tiers, solvers))
    if args.show_notes:
        chunks.append(build_notes(challenge, extended))

    output = "\n\n".join(chunks).strip() + "\n"
    if args.output:
        out_path = Path(args.output).resolve()
        out_path.write_text(output)
        print(f"Wrote table output: {out_path}")
    else:
        print(output, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
