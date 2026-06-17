#!/usr/bin/env python3
"""
build_combined_pdf.py -- consolidate per-module .tex docs into a single PDF.

Reads each standalone module .tex in package_module_docs/, strips the
preamble/title/maketitle, and emits one master document with each module as a
\\chapter. The original .tex files stay untouched (they remain individually
compilable for standalone use).

Usage:
    python3 docs/latex/build_combined_pdf.py
    # produces docs/latex/MFGraphs-package-docs.pdf
"""

import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE_DIR = HERE / "package_module_docs"
OUT_DIR = HERE
OUT_BASENAME = "MFGraphs-package-docs"

# Ordered chapters. (basename without .tex, chapter title, part_header_or_None)
CHAPTERS = [
    ("MFGraphs_loader",    "MFGraphs.wl -- Top-level loader",      "Loader"),
    ("primitives",         "primitives.wl",                         "Active package surface"),
    ("utilities",          "utilities.wl",                          None),
    ("scenarioTools",      "scenarioTools.wl",                      None),
    ("examples",           "examples.wl",                           None),
    ("unknownsTools",      "unknownsTools.wl",                      None),
    ("systemTools",        "systemTools.wl",                        None),
    ("solversTools",       "solversTools.wl",                       None),
    ("orchestrationTools", "orchestrationTools.wl",                 None),
    ("graphicsTools",      "graphicsTools.wl",                      None),
    ("Tawaf",              "Tawaf.wl (opt-in)",                     "Opt-in subpackages"),
    ("numericOracle",      "numericOracle.wl (opt-in)",             None),
]

MASTER_PREAMBLE = r"""\documentclass[11pt,oneside]{report}
\usepackage[margin=1in]{geometry}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{amsmath, amssymb}
\usepackage{array}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{xcolor}
\usepackage[colorlinks=true, linkcolor=blue!50!black, urlcolor=blue!50!black]{hyperref}
\title{MFGraphs Package Documentation\\
       \large Per-module long-form reference}
\author{MFGraphs maintainers}
\date{}
\begin{document}
\maketitle
\tableofcontents
\clearpage
"""

MASTER_POSTAMBLE = r"\end{document}" + "\n"


def extract_body(tex_path: Path) -> str:
    """Return the content between \\begin{document} and \\end{document},
    with \\maketitle removed."""
    src = tex_path.read_text()
    m = re.search(r"\\begin\{document\}(.*?)\\end\{document\}", src, re.DOTALL)
    if not m:
        raise RuntimeError(f"no document body in {tex_path}")
    body = m.group(1)
    # Drop \maketitle (the master provides its own title).
    body = re.sub(r"\\maketitle\s*", "", body)
    return body.strip()


def build_master_tex() -> str:
    parts = [MASTER_PREAMBLE]
    for basename, title, part_header in CHAPTERS:
        tex_path = MODULE_DIR / f"{basename}.tex"
        if not tex_path.exists():
            raise RuntimeError(f"missing: {tex_path}")
        if part_header:
            parts.append(f"\n\\part{{{part_header}}}\n")
        parts.append(f"\n\\chapter{{{title}}}\n")
        parts.append(extract_body(tex_path))
        parts.append("\n")
    parts.append(MASTER_POSTAMBLE)
    return "".join(parts)


def compile_pdf(master_tex: str) -> Path:
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        master = tmp_path / f"{OUT_BASENAME}.tex"
        master.write_text(master_tex)
        for _ in range(2):  # twice for TOC
            res = subprocess.run(
                ["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                 master.name],
                cwd=tmp_path,
                capture_output=True, text=True,
            )
            if res.returncode != 0:
                sys.stderr.write(res.stdout[-4000:])
                sys.stderr.write(res.stderr[-2000:])
                raise RuntimeError("pdflatex failed")
        pdf_src = tmp_path / f"{OUT_BASENAME}.pdf"
        pdf_dst = OUT_DIR / f"{OUT_BASENAME}.pdf"
        pdf_dst.write_bytes(pdf_src.read_bytes())
        return pdf_dst


def main() -> int:
    master_tex = build_master_tex()
    pdf_path = compile_pdf(master_tex)
    print(f"Wrote {pdf_path} ({pdf_path.stat().st_size:,} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
