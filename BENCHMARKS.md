# MFGraphs Performance Benchmarks

## Overview

This document describes the active benchmarking flow for the MFGraphs package.
The active benchmark targets the current scenario-kernel structural solver path:

```text
scenario -> makeSystem -> reduceSystem
```

Legacy tiered benchmarks for `DataToEquations` / `CriticalCongestionSolver`
are archived under `Scripts/archive/` and are not part of the default package
workflow while those solver-era symbols remain unloaded.

## Running benchmarks

### Prerequisites

- Wolfram Mathematica 12.0+ with `wolframscript` on PATH
- MFGraphs package (this repository)

### Active reduceSystem benchmark

Run the representative benchmark set:

```bash
wolframscript -file Scripts/BenchmarkReduceSystem.wls
```

Add a tag to append a dated history entry to this file:

```bash
wolframscript -file Scripts/BenchmarkReduceSystem.wls --tag "baseline"
```

Run a specific case with a custom timeout:

```bash
wolframscript -file Scripts/BenchmarkReduceSystem.wls --case example-12 --timeout 300
```

Results are exported to:

- `Results/reduce_system_latest.csv`
- `Results/reduce_system_<timestamp>.csv`
- `Results/reduce_system_solutions_latest.wl`
- `Results/reduce_system_solutions_<timestamp>.wl`

Current cases:

| Case | Description |
|------|-------------|
| `chain-2v` | two-vertex chain |
| `chain-3v-1exit` | three-vertex chain, one exit |
| `chain-3v-2exit` | three-vertex chain, two exits |
| `example-7` | built-in example 7 |
| `chain-5v-1exit` | five-vertex chain, one exit |
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
`BenchmarkReduceSystem.wls` workflow.

| Tier | Cases | Timeout | Description |
|------|-------|---------|-------------|
| small | 1-6, 27 | 60s | Linear chains, simple cycle |
| core | 9, 10, 12, 14, 15, 18, Paper | 300s | Stable representative cases for routine benchmarking |
| stress | 7, 8, 11, 17, 104, triangle-with-two-exits | 1800s | Symmetry-heavy or solver-stress cases; opt in explicitly |
| large | 13, 19-23, Braess variants, Jamaratv9, Grid0303 | 900s | Multi-entrance/exit, larger graphs |
| vlarge | Grid1020 | 1800s | 200-vertex grid; stress test designed to hit RecursionLimit |

## Solver benchmarked

**reduceSystem** is the current symbolic structural-system solver
(`MFGraphs/solversTools.wl`) over the reals. It supports critical congestion
systems only (`Alpha == 1` on every edge) and returns either rules for a fully
determined solution or a rules-plus-residual association for underdetermined
systems.

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
| `Scripts/BenchmarkReduceSystem.wls` | Active representative benchmark for `reduceSystem` |
| `Scripts/BenchmarkPreprocessing.wls` | Archived-preprocessing benchmark helper |
| `Scripts/perf_review_targeted.wls` | Targeted performance review helper |
| `Scripts/archive/BenchmarkSuite.wls` | Legacy tiered solver benchmark |
| `Scripts/archive/BottleneckReport.wls` | Legacy bottleneck profiler |

## Results format

Benchmark results are exported as CSV and JSON with these fields:

| Field | Description |
|-------|-------------|
| Case | Test case identifier |
| BuildMs | `makeSystem` construction time in milliseconds |
| WarmupMs | First `reduceSystem` run in milliseconds |
| RepMs | Repeated-timing estimate in milliseconds |
| Status | `OK` or `TIMEOUT` |
| Kind | `Rules`, `Underdetermined`, `NoSolution`, or `Other` |
| Valid | `isValidSystemSolution` result |

## Recommendations for future work

1. Broaden `BenchmarkReduceSystem.wls` with additional scenario-kernel cases once they are stable in the active API.
2. Keep legacy benchmark scripts archived until their dependencies are restored to the default package load path.
3. Track solver regressions through tagged benchmark entries plus active test-suite runs.
