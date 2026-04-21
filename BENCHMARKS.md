# MFGraphs Performance Benchmarks

## Overview

This document describes the active benchmarking flow for the MFGraphs package.
The benchmark suite now targets the critical solver path only.

## Running benchmarks

### Prerequisites

- Wolfram Mathematica 12.0+ with `wolframscript` on PATH
- MFGraphs package (this repository)

### Benchmark suite

Run the full benchmark suite:

```bash
wolframscript -file Scripts/BenchmarkSuite.wls
```

Run a specific tier only:

```bash
wolframscript -file Scripts/BenchmarkSuite.wls small
wolframscript -file Scripts/BenchmarkSuite.wls core
wolframscript -file Scripts/BenchmarkSuite.wls stress
wolframscript -file Scripts/BenchmarkSuite.wls large
wolframscript -file Scripts/BenchmarkSuite.wls vlarge
```

Run a specific case with a custom timeout:

```bash
wolframscript -file Scripts/BenchmarkSuite.wls core case=7 timeout=3600
```

Results are exported to `Results/benchmark_latest.csv` and `Results/benchmark_latest.json`, plus timestamped copies.

### Bottleneck profiling

Run the profiling analysis:

```bash
wolframscript -file Scripts/BottleneckReport.wls
```

This generates `Results/bottleneck_report.md` with detailed call counts and timing breakdowns for each profiled function.

## Test case tiers

| Tier | Cases | Timeout | Description |
|------|-------|---------|-------------|
| small | 1-6, 27 | 60s | Linear chains, simple cycle |
| core | 9, 10, 12, 14, 15, 18, Paper | 300s | Stable representative cases for routine benchmarking |
| stress | 7, 8, 11, 17, 104, triangle-with-two-exits | 1800s | Symmetry-heavy or solver-stress cases; opt in explicitly |
| large | 13, 19-23, Braess variants, Jamaratv9, Grid0303 | 900s | Multi-entrance/exit, larger graphs |
| vlarge | Grid1020 | 1800s | 200-vertex grid; stress test designed to hit RecursionLimit |

## Solver benchmarked

1. **CriticalCongestionSolver** -- Symbolic solver for the zero-flow case. Uses `Solve`, `Reduce`, `DNFReduce` (disjunctive normal form), and `TripleClean` (fixed-point simplification).

## Identified bottlenecks

### 1. DNFReduce exponential branching

`DNFReduce` in `DNFReduce.wl` recursively converts boolean expressions to disjunctive normal form. For `Or` expressions with N branches, this creates up to 2^N recursive paths. Each path calls `Solve` (via `ReDNFReduce`) and `Reduce`, with no caching of results across branches.

### 2. TripleClean fixed-point iteration

`TripleClean` in `DataToEquations.wl` repeatedly calls `TripleStep` until convergence. Each step calls `Solve` on the equality subsystem. Called multiple times per solver invocation (once in `MFGPreprocessing`, again in `MFGSystemSolver`).

## Optimizations implemented

### Optimization 1: DNFReduce Solve memoization

**File**: `MFGraphs/DNFReduce.wl`

Added a hash-based cache (`$SolveCache`) for `Solve` results within `ReDNFReduce`.

### Optimization 2: DNFReduce branch pruning

**File**: `MFGraphs/DNFReduce.wl`

Added early-exit checks and reduced branching work for trivial false branches.

### Optimization 3: TripleStep Solve memoization

**File**: `MFGraphs/DataToEquations.wl`

`TripleStep` now uses cached solve results for repeated equality systems.

## Profiling scripts

| Script | Purpose |
|--------|---------|
| `Scripts/BenchmarkSuite.wls` | Full benchmark across all tiers for `CriticalCongestion` |
| `Scripts/BenchmarkHelpers.wls` | Helper functions (SafeExecute, GraphMetadata, parameter tables) |
| `Scripts/ProfileInstrument.wls` | Non-invasive DownValues-based profiling wrappers |
| `Scripts/BottleneckReport.wls` | Runs profiled solver on representative cases |

## Results format

Benchmark results are exported as CSV and JSON with these fields:

| Field | Description |
|-------|-------------|
| Key | Test case identifier |
| Tier | small/core/stress/large/vlarge |
| Solver | CriticalCongestion |
| NumVertices | Graph vertex count |
| NumEdges | Graph edge count |
| D2ETime | DataToEquations wall time (seconds) |
| SolveTime | Solver wall time (seconds) |
| SolveMemory | Peak memory delta (bytes) |
| Status | OK/TIMEOUT/FAILED/SKIPPED |
| Residual | Infinity-norm of the shared Kirchhoff residual |
| EquationResidual | Infinity-norm of the equation residual, computed on the public solver result |
| StopReason | Convergence stop reason when reported |
| ConvergenceResidual | Reported convergence residual when present |
| Iterations | Iteration or accepted-step count reported by the solver |

## Recommendations for future work

1. **Parallel DNFReduce branches**: `Or` branches in `DNFReduce` are independent and could be evaluated in parallel.
2. **Sparse matrix operations**: For large graphs, prioritize sparse `LinearSolve` paths when possible.
