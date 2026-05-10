# MFGraphs Performance Benchmarks

## Overview

This document describes the active benchmarking flow for the MFGraphs package.
The active benchmark targets the current scenario-kernel structural solver path:

```text
scenario -> makeSystem -> dnfReduceSystem
```

Opt-in exact solver comparisons are available with `--solver optimizeddnf`
and `--solver activeset`; the default benchmark path remains `dnfReduceSystem`.

Legacy tiered benchmarks for `DataToEquations` / `CriticalCongestionSolver`
are archived under `Scripts/archive/` and are not part of the default package
workflow while those solver-era symbols remain unloaded.

## Benchmark policy for solver-sensitive changes

Any pull request that modifies solver logic, system construction, or preprocessing
must include a before/after benchmark run using `--tag`:

```bash
# Before the change (on the base branch)
wolframscript -file Scripts/BenchmarkSystemSolver.wls --tag "before: <description>" --timeout 60

# After the change (on the feature branch)
wolframscript -file Scripts/BenchmarkSystemSolver.wls --tag "after: <description>" --timeout 60
```

Include the resulting table rows from `BENCHMARKS.md` in the PR description.
Regressions above ~25% on non-timeout cases require explanation before merge.

Solver-sensitive pull requests must also include `ProfileScenarioKernel.wls` output
for at least one easy case (e.g. `chain-3v-1exit`) and one hard case (e.g. `grid-3x2`
or `example-12`). Record which phase absorbed the cost change: scenario build, unknowns,
system construction, solver inputs, or solve.

## Running benchmarks

### Prerequisites

- Wolfram Mathematica 12.0+ with `wolframscript` on PATH
- MFGraphs package (this repository)

### Active system solver benchmark

Run the representative benchmark set:

```bash
wolframscript -file Scripts/BenchmarkSystemSolver.wls
```

Add a tag to append a dated history entry to this file:

```bash
wolframscript -file Scripts/BenchmarkSystemSolver.wls --tag "baseline"
```

Run a specific case with a custom timeout:

```bash
wolframscript -file Scripts/BenchmarkSystemSolver.wls --case example-12 --timeout 60
```

Compare an alternate solver explicitly:

```bash
wolframscript -file Scripts/BenchmarkSystemSolver.wls --solver reduce --case example-12 --timeout 60
wolframscript -file Scripts/BenchmarkSystemSolver.wls --solver optimizeddnf --case example-12 --timeout 60
wolframscript -file Scripts/BenchmarkSystemSolver.wls --solver activeset --case example-12 --timeout 60
```

Results are exported to:

- `Results/system_solver_<solver>_latest.csv`
- `Results/system_solver_<solver>_<timestamp>.csv`
- `Results/system_solver_<solver>_solutions_latest.wl`
- `Results/system_solver_<solver>_solutions_<timestamp>.wl`

For non-mutating staged diagnostics, run:

```bash
wolframscript -file Scripts/ProfileScenarioKernel.wls --case example-12 --timeout 10
wolframscript -file Scripts/ProfileScenarioKernel.wls --solver optimizeddnf --case example-12 --timeout 10
```

Current cases:

| Case | Description |
|------|-------------|
| `chain-2v` | two-vertex chain |
| `chain-3v-1exit` | three-vertex chain, one exit |
| `chain-3v-2exit` | three-vertex chain, two exits |
| `example-7` | built-in example 7 |
| `chain-5v-1exit` | five-vertex chain, one exit |
| `grid-2x3` | 2 by 3 grid, entry `{1, 100}`, exit `{6, 0}` |
| `grid-3x2` | 3 by 2 grid, entry `{1, 100}`, exit `{6, 0}` |
| `example-12` | built-in example 12 |

### Archived solver-era benchmarks

The old tiered benchmark runner and bottleneck profiler are retained only for
manual legacy work:

```bash
wolframscript -file Scripts/archive/BenchmarkSuite.wls
wolframscript -file Scripts/archive/BottleneckReport.wls
```

See `Scripts/LEGACY_SOLVER_SCRIPTS.md` and
`Scripts/LEGACY_GETEXAMPLEDATA.md` before using them.

## Legacy test case tiers

These tiers describe `Scripts/archive/BenchmarkSuite.wls`, not the active
`BenchmarkSystemSolver.wls` workflow.

| Tier | Cases | Timeout | Description |
|------|-------|---------|-------------|
| small | 1-6, 27 | 60s | Linear chains, simple cycle |
| core | 9, 10, 12, 14, 15, 18, Paper | 300s | Stable representative cases for routine benchmarking |
| stress | 7, 8, 11, 17, 104, triangle-with-two-exits | 1800s | Symmetry-heavy or solver-stress cases; opt in explicitly |
| large | 13, 19-23, Braess variants, Jamaratv9, Grid0303 | 900s | Multi-entrance/exit, larger graphs |
| vlarge | Grid1020 | 1800s | 200-vertex grid; stress test designed to hit RecursionLimit |

## Solver benchmarked

**dnfReduceSystem** is the default symbolic structural-system solver
(`MFGraphs/solversTools.wl`). It uses the shared linear preprocessing path and
then applies equality-substitution plus disjunction distribution to avoid raw
`Reduce` timeouts. Raw **reduceSystem** remains available as a direct Wolfram
`Reduce` baseline.

**optimizedDNFReduceSystem** and **activeSetReduceSystem** are opt-in exact
solvers for critical-congestion systems. They carry small branch families
directly as rules plus residual constraints, fall back to the proven exact DNF
path for larger residual variable sets, and return the same rule/residual shape
as `dnfReduceSystem`.

## system solver benchmark history

DNF-first benchmark entries from `Scripts/BenchmarkSystemSolver.wls` are appended
here when the script is run with `--tag`.

### 2026-05-10 - post-cache

**Commit:** `bf53e8f`  
**Environment:** local `wolframscript`, `$MFGraphsVerbose=False`  
**Solver:** `dnf`  
**Timeout per case:** 60s  
**Solutions file:** `Results/system_solver_dnf_solutions_20260510-124613.wl`

| Case | Solver | Build (ms) | Warmup (ms) | Repeat (ms) | Status | Kind | Transition flows | Residual Jts | Valid |
|------|--------|-----------|------------|--------------|--------|------|------------------|--------------|-------|
| chain-2v | dnf | 8.45 | 1.03 | 0.74 | OK | Rules | Unique | 0 | True |
| chain-3v-1exit | dnf | 2.42 | 1.01 | 0.99 | OK | Rules | Unique | 0 | True |
| chain-3v-2exit | dnf | 2.32 | 29.72 | 2.85 | OK | Underdetermined | Unique | 0 | True |
| chain-3-midentry | dnf | 3.08 | 3.77 | 3.29 | OK | Rules | Unique | 0 | True |
| chain-5v-1exit | dnf | 3.25 | 1.64 | 1.51 | OK | Rules | Unique | 0 | True |
| grid-2x3 | dnf | 10.93 | 7.53 | 7.08 | OK | Rules | Unique | 0 | True |
| grid-3x2 | dnf | 10.33 | 7.54 | 6.86 | OK | Rules | Unique | 0 | True |
| example-12 | dnf | 7.22 | 8.43 | 8.2 | OK | Rules | Unique | 0 | True |
---

### 2026-05-10 - post-rowreduce-opt

**Commit:** `bf53e8f`  
**Environment:** local `wolframscript`, `$MFGraphsVerbose=False`  
**Solver:** `optimizeddnf`  
**Timeout per case:** 60s  
**Solutions file:** `Results/system_solver_optimizeddnf_solutions_20260510-123459.wl`

| Case | Solver | Build (ms) | Warmup (ms) | Repeat (ms) | Status | Kind | Transition flows | Residual Jts | Valid |
|------|--------|-----------|------------|--------------|--------|------|------------------|--------------|-------|
| chain-2v | optimizeddnf | 8.53 | 1.19 | 0.97 | OK | Rules | Unique | 0 | True |
| chain-3v-1exit | optimizeddnf | 2.63 | 1.1 | 1.15 | OK | Rules | Unique | 0 | True |
| chain-3v-2exit | optimizeddnf | 2.87 | 30.29 | 4.75 | OK | Underdetermined | Underdetermined | 1 | True |
| chain-3-midentry | optimizeddnf | 3.15 | 7.69 | 5.8 | OK | Underdetermined | Underdetermined | 1 | True |
| chain-5v-1exit | optimizeddnf | 3.77 | 1.68 | 1.54 | OK | Rules | Unique | 0 | True |
| grid-2x3 | optimizeddnf | 11.81 | 9.57 | 7.8 | OK | Rules | Unique | 0 | True |
| grid-3x2 | optimizeddnf | 12.05 | 7.54 | 7.75 | OK | Rules | Unique | 0 | True |
| example-12 | optimizeddnf | 8.19 | 10.22 | 9.28 | OK | Rules | Unique | 0 | True |
---

### 2026-05-10 - post-rowreduce

**Commit:** `bf53e8f`  
**Environment:** local `wolframscript`, `$MFGraphsVerbose=False`  
**Solver:** `dnf`  
**Timeout per case:** 60s  
**Solutions file:** `Results/system_solver_dnf_solutions_20260510-123456.wl`

| Case | Solver | Build (ms) | Warmup (ms) | Repeat (ms) | Status | Kind | Transition flows | Residual Jts | Valid |
|------|--------|-----------|------------|--------------|--------|------|------------------|--------------|-------|
| chain-2v | dnf | 8.51 | 1.06 | 0.77 | OK | Rules | Unique | 0 | True |
| chain-3v-1exit | dnf | 2.49 | 1.16 | 1.15 | OK | Rules | Unique | 0 | True |
| chain-3v-2exit | dnf | 2.83 | 30.65 | 2.92 | OK | Underdetermined | Unique | 0 | True |
| chain-3-midentry | dnf | 3.06 | 3.76 | 3.37 | OK | Rules | Unique | 0 | True |
| chain-5v-1exit | dnf | 3.54 | 1.58 | 1.66 | OK | Rules | Unique | 0 | True |
| grid-2x3 | dnf | 11.43 | 7.52 | 7.51 | OK | Rules | Unique | 0 | True |
| grid-3x2 | dnf | 10.42 | 7.13 | 6.88 | OK | Rules | Unique | 0 | True |
| example-12 | dnf | 7.2 | 8.27 | 8.07 | OK | Rules | Unique | 0 | True |
---

### 2026-05-10 - post-dnf

**Commit:** `bf53e8f`  
**Environment:** local `wolframscript`, `$MFGraphsVerbose=False`  
**Solver:** `dnf`  
**Timeout per case:** 60s  
**Solutions file:** `Results/system_solver_dnf_solutions_20260510-115527.wl`

| Case | Solver | Build (ms) | Warmup (ms) | Repeat (ms) | Status | Kind | Transition flows | Residual Jts | Valid |
|------|--------|-----------|------------|--------------|--------|------|------------------|--------------|-------|
| chain-2v | dnf | 8.24 | 23.58 | 0.72 | OK | Rules | Unique | 0 | True |
| chain-3v-1exit | dnf | 2.7 | 0.89 | 0.81 | OK | Rules | Unique | 0 | True |
| chain-3v-2exit | dnf | 2.42 | 7.04 | 2.74 | OK | Underdetermined | Unique | 0 | True |
| chain-3-midentry | dnf | 3.02 | 3.68 | 3.24 | OK | Rules | Unique | 0 | True |
| chain-5v-1exit | dnf | 3.61 | 1.38 | 1.31 | OK | Rules | Unique | 0 | True |
| grid-2x3 | dnf | 11.2 | 7.3 | 7.1 | OK | Rules | Unique | 0 | True |
| grid-3x2 | dnf | 10.9 | 7. | 7.33 | OK | Rules | Unique | 0 | True |
| example-12 | dnf | 7.28 | 8.7 | 8.87 | OK | Rules | Unique | 0 | True |
---

### 2026-05-10 - pre-dnf

**Commit:** `bf53e8f`  
**Environment:** local `wolframscript`, `$MFGraphsVerbose=False`  
**Solver:** `dnf`  
**Timeout per case:** 60s  
**Solutions file:** `Results/system_solver_dnf_solutions_20260510-115524.wl`

| Case | Solver | Build (ms) | Warmup (ms) | Repeat (ms) | Status | Kind | Transition flows | Residual Jts | Valid |
|------|--------|-----------|------------|--------------|--------|------|------------------|--------------|-------|
| chain-2v | dnf | 8.65 | 24.15 | 0.73 | OK | Rules | Unique | 0 | True |
| chain-3v-1exit | dnf | 2.8 | 0.94 | 0.96 | OK | Rules | Unique | 0 | True |
| chain-3v-2exit | dnf | 2.35 | 7.6 | 2.8 | OK | Underdetermined | Unique | 0 | True |
| chain-3-midentry | dnf | 3.1 | 3.69 | 3. | OK | Rules | Unique | 0 | True |
| chain-5v-1exit | dnf | 3.4 | 1.37 | 1.36 | OK | Rules | Unique | 0 | True |
| grid-2x3 | dnf | 11.45 | 7.28 | 6.92 | OK | Rules | Unique | 0 | True |
| grid-3x2 | dnf | 10.64 | 7.51 | 6.88 | OK | Rules | Unique | 0 | True |
| example-12 | dnf | 7.41 | 9.06 | 8.47 | OK | Rules | Unique | 0 | True |
---

### 2026-05-10 - post-rename

**Commit:** `bf53e8f`  
**Environment:** local `wolframscript`, `$MFGraphsVerbose=False`  
**Solver:** `optimizeddnf`  
**Timeout per case:** 60s  
**Solutions file:** `Results/system_solver_optimizeddnf_solutions_20260510-115448.wl`

| Case | Solver | Build (ms) | Warmup (ms) | Repeat (ms) | Status | Kind | Transition flows | Residual Jts | Valid |
|------|--------|-----------|------------|--------------|--------|------|------------------|--------------|-------|
| chain-2v | optimizeddnf | 8.36 | 23.84 | 0.68 | OK | Rules | Unique | 0 | True |
| chain-3v-1exit | optimizeddnf | 2.52 | 0.86 | 0.86 | OK | Rules | Unique | 0 | True |
| chain-3v-2exit | optimizeddnf | 2.29 | 26.97 | 4.34 | OK | Underdetermined | Underdetermined | 1 | True |
| chain-3-midentry | optimizeddnf | 3.25 | 7.83 | 5.52 | OK | Underdetermined | Underdetermined | 1 | True |
| chain-5v-1exit | optimizeddnf | 3.72 | 1.44 | 1.31 | OK | Rules | Unique | 0 | True |
| grid-2x3 | optimizeddnf | 11.35 | 7.25 | 7.35 | OK | Rules | Unique | 0 | True |
| grid-3x2 | optimizeddnf | 11.18 | 7.42 | 6.67 | OK | Rules | Unique | 0 | True |
| example-12 | optimizeddnf | 7.11 | 8.75 | 8.72 | OK | Rules | Unique | 0 | True |
---

### 2026-05-10 - pre-rename

**Commit:** `bf53e8f`  
**Environment:** local `wolframscript`, `$MFGraphsVerbose=False`  
**Solver:** `optimizeddnf`  
**Timeout per case:** 60s  
**Solutions file:** `Results/system_solver_optimizeddnf_solutions_20260510-115422.wl`

| Case | Solver | Build (ms) | Warmup (ms) | Repeat (ms) | Status | Kind | Transition flows | Residual Jts | Valid |
|------|--------|-----------|------------|--------------|--------|------|------------------|--------------|-------|
| chain-2v | optimizeddnf | 8.57 | 23.94 | 0.7 | OK | Rules | Unique | 0 | True |
| chain-3v-1exit | optimizeddnf | 2.62 | 0.86 | 0.89 | OK | Rules | Unique | 0 | True |
| chain-3v-2exit | optimizeddnf | 2.49 | 26.64 | 3.92 | OK | Underdetermined | Underdetermined | 1 | True |
| chain-3-midentry | optimizeddnf | 3.26 | 6.92 | 4.7 | OK | Underdetermined | Underdetermined | 1 | True |
| chain-5v-1exit | optimizeddnf | 3.61 | 1.48 | 1.34 | OK | Rules | Unique | 0 | True |
| grid-2x3 | optimizeddnf | 11.33 | 7.58 | 7.12 | OK | Rules | Unique | 0 | True |
| grid-3x2 | optimizeddnf | 11.1 | 6.81 | 6.46 | OK | Rules | Unique | 0 | True |
| example-12 | optimizeddnf | 7.11 | 8.42 | 8.38 | OK | Rules | Unique | 0 | True |
---

### 2026-05-07 - post-PR-192

**Commit:** `f0a4530`  
**Environment:** local `wolframscript`, `$MFGraphsVerbose=False`  
**Solver:** `dnf`  
**Timeout per case:** 60s  
**Solutions file:** `Results/system_solver_dnf_solutions_20260507-005846.wl`

| Case | Solver | Build (ms) | Warmup (ms) | Repeat (ms) | Status | Kind | Transition flows | Residual Jts | Valid |
|------|--------|-----------|------------|--------------|--------|------|------------------|--------------|-------|
| chain-2v | dnf | 9.8 | 24.73 | 0.67 | OK | Rules | Unique | 0 | True |
| chain-3v-1exit | dnf | 2.57 | 0.93 | 0.88 | OK | Rules | Unique | 0 | True |
| chain-3v-2exit | dnf | 2.31 | 7.54 | 2.69 | OK | Underdetermined | Unique | 0 | True |
| example-7 | dnf | 3.63 | 4.82 | 4.67 | OK | Rules | Unique | 0 | True |
| chain-5v-1exit | dnf | 4.47 | 1.89 | 1.6 | OK | Rules | Unique | 0 | True |
| grid-2x3 | dnf | 12.54 | 7.79 | 7.01 | OK | Rules | Unique | 0 | True |
| grid-3x2 | dnf | 11.27 | 7.07 | 6.74 | OK | Rules | Unique | 0 | True |
| example-12 | dnf | 5.09 | 8.51 | 8.38 | OK | Rules | Unique | 0 | True |
---

### 2026-05-05 - unconditional EqGeneral

**Commit:** `4648fac`  
**Environment:** local `wolframscript`, `$MFGraphsVerbose=False`  
**Solver:** `dnf`  
**Timeout per case:** 60s  
**Solutions file:** `Results/system_solver_dnf_solutions_20260505-112802.wl`

| Case | Solver | Build (ms) | Warmup (ms) | Repeat (ms) | Status | Kind | Transition flows | Residual Jts | Valid |
|------|--------|-----------|------------|--------------|--------|------|------------------|--------------|-------|
| chain-2v | dnf | 9.3 | 25.6 | 0.78 | OK | Rules | Unique | 0 | True |
| chain-3v-1exit | dnf | 3.11 | 1.37 | 1.07 | OK | Rules | Unique | 0 | True |
| chain-3v-2exit | dnf | 2.7 | 2.1 | 2.02 | OK | Rules | Unique | 0 | True |
| example-7 | dnf | 4.36 | 2.79 | 2.62 | OK | Rules | Unique | 0 | True |
| chain-5v-1exit | dnf | 3.49 | 1.8 | 1.64 | OK | Rules | Unique | 0 | True |
| grid-2x3 | dnf | 11.12 | 6.76 | 6.43 | OK | Rules | Unique | 0 | True |
| grid-3x2 | dnf | 10.91 | 6.98 | 6.28 | OK | Rules | Unique | 0 | True |
| example-12 | dnf | 5.11 | 7.62 | 7.53 | OK | Rules | Unique | 0 | True |
---

### 2026-05-05 - post-ship-20260505

**Commit:** `4648fac`  
**Environment:** local `wolframscript`, `$MFGraphsVerbose=False`  
**Solver:** `dnf`  
**Timeout per case:** 60s  
**Solutions file:** `Results/system_solver_dnf_solutions_20260505-102048.wl`

| Case | Solver | Build (ms) | Warmup (ms) | Repeat (ms) | Status | Kind | Transition flows | Residual Jts | Valid |
|------|--------|-----------|------------|--------------|--------|------|------------------|--------------|-------|
| chain-2v | dnf | 8.63 | 25.39 | 0.71 | OK | Rules | Unique | 0 | True |
| chain-3v-1exit | dnf | 2.56 | 0.98 | 0.9 | OK | Rules | Unique | 0 | True |
| chain-3v-2exit | dnf | 2.28 | 2.33 | 1.95 | OK | Rules | Unique | 0 | True |
| example-7 | dnf | 4.18 | 4.38 | 3.96 | OK | Rules | Unique | 0 | True |
| chain-5v-1exit | dnf | 3.41 | 1.57 | 1.41 | OK | Rules | Unique | 0 | True |
| grid-2x3 | dnf | 10.89 | 425.63 | 438.73 | OK | Branched | Underdetermined | 4 | True |
| grid-3x2 | dnf | 13. | 437.22 | 425.26 | OK | Branched | Underdetermined | 4 | True |
| example-12 | dnf | 6.24 | 268.24 | 213.85 | OK | Branched | Underdetermined | 3 | True |
---

### 2026-05-04 - issue6-list-storage

**Commit:** `045db36`  
**Environment:** local `wolframscript`, `$MFGraphsVerbose=False`  
**Solver:** `dnf`  
**Timeout per case:** 60s  
**Solutions file:** `Results/system_solver_dnf_solutions_20260504-102432.wl`

| Case | Solver | Build (ms) | Warmup (ms) | Repeat (ms) | Status | Kind | Transition flows | Residual Jts | Valid |
|------|--------|-----------|------------|--------------|--------|------|------------------|--------------|-------|
| chain-2v | dnf | 9.47 | 25.31 | 0.68 | OK | Rules | Unique | 0 | True |
| chain-3v-1exit | dnf | 2.5 | 0.94 | 0.84 | OK | Rules | Unique | 0 | True |
| chain-3v-2exit | dnf | 2.49 | 2.59 | 1.97 | OK | Rules | Unique | 0 | True |
| example-7 | dnf | 4.08 | 4.32 | 4.02 | OK | Rules | Unique | 0 | True |
| chain-5v-1exit | dnf | 3.75 | 1.62 | 1.45 | OK | Rules | Unique | 0 | True |
| grid-2x3 | dnf | 11. | 628.09 | 576.54 | OK | Branched | Underdetermined | 1 | True |
| grid-3x2 | dnf | 10.84 | 548.82 | 524.39 | OK | Branched | Unique | 0 | True |
| example-12 | dnf | 5.37 | 312.77 | 285.24 | OK | Branched | Underdetermined | 1 | True |
---

### 2026-04-30 — post-PR-168 baseline (dnf, 60s timeout)

**Commit:** `0edaaa6` (Merge PR #168)
**Environment:** `wolframscript`, Mathematica 14.3.0 ARM (Mac), `$MFGraphsVerbose=False`

| Case | Build (ms) | Warmup (ms) | Repeat (ms) | Kind | Valid |
|------|-----------|------------|------------|------|-------|
| chain-2v | 9.3 | 25.8 | 0.95 | Rules | True |
| chain-3v-1exit | 2.81 | 1.19 | 1.09 | Rules | True |
| chain-3v-2exit | 2.32 | 31.79 | 3.11 | Branched | True |
| example-7 | 3.87 | 13.23 | 7.02 | Branched | True |
| chain-5v-1exit | 3.42 | 2.02 | 1.76 | Rules | True |
| grid-2x3 | 10.95 | 12846 | **11417** | Branched | True |
| grid-3x2 | 10.82 | 3096 | **2577** | Branched | True |
| example-12 | 5.23 | 1797 | **1543** | Branched | True |

Notes:
- `grid-2x3` and `grid-3x2` have identical vertex/edge counts but the DNF solve time differs by ~4.4×, reflecting topology sensitivity in the disjunction ordering. The column-major layout (2×3) produces significantly more branch combinations than the row-major layout (3×2).
- `example-12` settles at ~1.5s repeat; `reduceSystem` times out at 160s on this case (see history below).

## reduceSystem benchmark history (manual)

### 2026-04-25 — core notebook solver sanity + timeout sweep

**Commit:** `e7ed719`  
**Environment:** local `wolframscript`, `$MFGraphsVerbose=False`

Method:
- Build scenario -> `makeSystem` -> `reduceSystem`.
- For baseline cases: one warmup + `RepeatedTiming[..., 5]`.
- For hard case (Example 12): `TimeConstrained` sweep at 20/40/80/160 seconds.

#### Baseline cases

| Case | Warmup (ms) | Repeated (ms) | Head | NonFalse |
|------|-------------|---------------|------|----------|
| 2-vertex chain | 317.07 | 2.04 | `And` | True |
| `GridScenario[{3}]` one-exit | 5.84 | 5.43 | `And` | True |
| `GridScenario[{3}]` two-exit | 13.09 | 12.74 | `And` | True |
| `Example 7` | 158.34 | 157.81 | `And` | True |
| `GridScenario[{5}]` one-exit | 284.71 | 283.87 | `And` | True |

#### Timeout sweep — `Example 12`

| Timeout (s) | Elapsed (s) | TimedOut | Head |
|-------------|-------------|----------|------|
| 20 | 20.13 | True | `TimedOut` |
| 40 | 40.13 | True | `TimedOut` |
| 80 | 82.33 | True | `TimedOut` |
| 160 | 168.68 | True | `TimedOut` |

Interpretation:
- `reduceSystem` is fast on small/core examples and moderate on larger chains.
- `Example 12` does not complete even at 160s under current formulation.

## Related scripts

| Script | Purpose |
|--------|---------|
| `Scripts/BenchmarkSystemSolver.wls` | Active DNF-first benchmark; supports `--solver` for comparisons |
| `Scripts/ProfileScenarioKernel.wls` | Non-mutating staged profiler for construction, preprocessing, and solver time |
| `Scripts/ProfileDNFReduce.wls` | DNF reducer diagnostic profiler; sweeps ordering strategies and reports branch/substitution metrics |
| `Scripts/BenchmarkReduceSystem.wls` | Historical raw-Reduce baseline benchmark |
| `Scripts/BenchmarkPreprocessing.wls` | Archived-preprocessing benchmark helper |
| `Scripts/perf_review_targeted.wls` | Targeted performance review helper |
| `Scripts/archive/BenchmarkSuite.wls` | Legacy tiered solver benchmark |
| `Scripts/archive/BottleneckReport.wls` | Legacy bottleneck profiler |

## Results format

Benchmark results are exported as CSV and JSON with these fields:

| Field | Description |
|-------|-------------|
| Case | Test case identifier |
| Solver | Solver name (`dnf`, `optimizeddnf`, `activeset`, `reduce`, `boolean`, or `findinstance`) when using `BenchmarkSystemSolver.wls` |
| BuildMs | `makeSystem` construction time in milliseconds |
| WarmupMs | First solver run in milliseconds |
| RepMs | One timed repeat after warmup, in milliseconds |
| Status | `OK` or `TIMEOUT` |
| Kind | `Rules` (fully determined), `Branched` (residual with top-level disjunctions), `Parametric` (tracked residual without branching), `NoSolution`, or `Unknown` — as classified by `solutionResultKind` |
| Valid | `isValidSystemSolution` result |

## Recommendations for future work

1. Broaden `BenchmarkSystemSolver.wls` with additional scenario-kernel cases once they are stable in the active API.
2. Keep legacy benchmark scripts archived until their dependencies are restored to the default package load path.
3. Track solver regressions through tagged benchmark entries plus active test-suite runs.