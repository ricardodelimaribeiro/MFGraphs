# MFGraphs Performance Benchmarks

## Overview

This document describes the benchmarking and profiling infrastructure for the MFGraphs package, along with the algorithmic optimizations implemented to address identified bottlenecks.

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
wolframscript -file Scripts/BenchmarkSuite.wls medium
wolframscript -file Scripts/BenchmarkSuite.wls large
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
| medium | 7-15, 17-18, 104, triangle, Paper | 300s | Y-networks, triangles, switching costs |
| large | 13, 19-23, Braess variants, Jamaratv9, Grid0303 | 900s | Multi-entrance/exit, larger graphs |
| vlarge | Grid1020 | 1800s | 200-vertex grid (symbolic solvers may not terminate) |

## Solvers benchmarked

1. **CriticalCongestionSolver** -- Symbolic solver for the zero-flow case. Uses `Solve`, `Reduce`, `DNFReduce` (disjunctive normal form), and `TripleClean` (fixed-point simplification).

2. **NonLinearSolver** -- Iterative fixed-point solver for general congestion. Calls `MFGSystemSolver` up to 15 times per solve, with `FindRoot`/`NIntegrate` for Hamiltonian cost evaluation.

3. **MonotoneSolver** -- ODE-based gradient flow using `NDSolve` with `GradientProjection` (PseudoInverse-based).

## Identified bottlenecks

### 1. DNFReduce exponential branching

`DNFReduce` in `DNFReduce.m` recursively converts boolean expressions to disjunctive normal form. For `Or` expressions with N branches, this creates up to 2^N recursive paths. Each path calls `Solve` (via `ReDNFReduce`) and `Reduce`, with no caching of results across branches.

**Impact**: Dominates CriticalCongestionSolver runtime for cases with switching costs (keys 8, 10, 11, 14, Braess variants).

### 2. TripleClean fixed-point iteration

`TripleClean` in `DataToEquations.m` repeatedly calls `TripleStep` until convergence. Each step calls `Solve` on the equality subsystem. Called multiple times per solver invocation (once in `MFGPreprocessing`, again in `MFGSystemSolver`).

**Impact**: 5-10 iterations typical; each iteration involves symbolic solve of increasing complexity.

### 3. PseudoInverse in GradientProjection

`GradientProjection` in `Monotone.m` computes `PseudoInverse` (O(n^3)) at every ODE evaluation step. For `NDSolve` with 100 time steps and adaptive stepping, this can mean hundreds of PseudoInverse calls with nearly-identical matrices.

**Impact**: Dominates MonotoneSolver runtime for larger graphs (Grid0303+).

### 4. Per-point FindRoot/NIntegrate

`M[j, x, edge]` in `NonLinearSolver.m` calls `FindRoot` for every (j, x) evaluation point. `IntegratedMass` wraps this in `NIntegrate`, which samples M at many points. Called per-edge, per-iteration of `NonLinearSolver`.

**Impact**: Significant for problems with many edges and many `NonLinearSolver` iterations.

## Optimizations implemented

### Optimization 1: DNFReduce Solve memoization

**File**: `MFGraphs/DNFReduce.m`

Added a hash-based cache (`$SolveCache`) for `Solve` results within `ReDNFReduce`. The same equality often appears across multiple branches of an `Or` expression; memoization avoids redundant symbolic solves.

- `CachedSolve[eq]` -- hash-based lookup/store wrapper around `Solve`
- `ClearSolveCache[]` -- called at solver entry points to prevent stale results

**Expected impact**: 30-70% fewer `Solve` calls on cases with switching costs.

### Optimization 2: DNFReduce branch pruning

**File**: `MFGraphs/DNFReduce.m`

Added early-exit checks:
- `DNFReduce[xp, And[fst, rst]]` returns `False` immediately when `xp === False`
- `DNFReduce[xp, Or[fst, scd]]` skips `False` branches in the result instead of creating trivial Or expressions

**Expected impact**: Avoids exponential blowup on inconsistent branches.

### Optimization 3: PseudoInverse caching in MonotoneSolver

**File**: `MFGraphs/Monotone.m`

Added `CachedGradientProjection` that caches the PseudoInverse matrix and reuses it when the flow vector `x` hasn't changed beyond a tolerance threshold (default 10^-6).

- Controlled by `"UseCachedGradient" -> True` option (default)
- Set `"UseCachedGradient" -> False` to use the original uncached behavior

**Expected impact**: Significant speedup for larger graphs where PseudoInverse dominates.

### Optimization 4: M interpolation for NonLinearSolver

**File**: `MFGraphs/NonLinearSolver.m`

Added `PrecomputeM[jMin, jMax, edge, nPoints]` that builds an `InterpolatingFunction` from a grid of `FindRoot` evaluations. `FastIntegratedMass` and `FastIntegratedMass` use this interpolation instead of per-point `FindRoot` calls.

- `nPoints` controls grid resolution (default 50)
- Users can switch between exact (`Cost`/`IntegratedMass`) and fast (`FastIntegratedMass`/`FastIntegratedMass`) versions

**Expected impact**: Large speedup for NIntegrate-heavy iterations; slight numerical approximation controlled by grid resolution.

### Optimization 5: DNFReduce sequential And-Or with early exit

**File**: `MFGraphs/DNFReduce.m`

The And-Or distribution case (`DNFReduce[xp_, And[fst_Or, rst_]]`) previously evaluated all disjuncts of `fst` unconditionally using `Map`. Replaced with a sequential `Do` loop and `Catch/Throw` that exits as soon as any branch returns `r === xp` — the same short-circuit condition used in the existing 2-arg Or case. When a branch leaves `xp` unchanged, it means the Or constraint is already satisfied, so remaining branches cannot contribute new information.

**Expected impact**: Fewer branch evaluations for And-Or patterns in well-constrained systems; symmetric to the existing Or short-circuit which showed dramatic speedup on cases 11 and 12.

### Optimization 6: TripleStep Solve memoization

**File**: `MFGraphs/DataToEquations.m`

Both overloads of `TripleStep` called `Solve` directly. Replaced with `CachedSolve` (already defined in `DNFReduce.m`) so that identical equality systems encountered within a single `MFGSystemSolver` invocation (between `ClearSolveCache` calls) are solved only once.

**Expected impact**: Modest reduction in `Solve` calls during multi-step `TripleClean` fixed-point iteration.

### Optimization 7: NonLinearSolver tolerance-based early stopping

**File**: `MFGraphs/NonLinearSolver.m`

Added `"Tolerance" -> 0` option to `NonLinearSolver`. When `Tolerance > 0`, iteration switches from `FixedPointList` (exact equality on Associations, rarely fires for floating-point results) to `NestWhileList` with a norm-based stopping condition: stops when the infinity-norm change in flow variables between consecutive iterations is below the tolerance.

**Usage**: `NonLinearSolver[d2e, "Tolerance" -> 10^-6]`

**Expected impact**: Significant savings when problems converge in fewer than `MaxIterations` steps; previously all 15 iterations ran regardless of numerical convergence.

## Profiling scripts

| Script | Purpose |
|--------|---------|
| `Scripts/BenchmarkSuite.wls` | Full benchmark across all tiers and solvers |
| `Scripts/BenchmarkHelpers.wls` | Helper functions (SafeExecute, GraphMetadata, parameter tables) |
| `Scripts/ProfileInstrument.wls` | Non-invasive DownValues-based profiling wrappers |
| `Scripts/BottleneckReport.wls` | Runs profiled solvers on representative cases |

## Results format

Benchmark results are exported as CSV and JSON with these fields:

| Field | Description |
|-------|-------------|
| Key | Test case identifier |
| Tier | small/medium/large/vlarge |
| Solver | CriticalCongestion/NonLinearSolver/Monotone |
| NumVertices | Graph vertex count |
| NumEdges | Graph edge count |
| D2ETime | DataToEquations wall time (seconds) |
| SolveTime | Solver wall time (seconds) |
| SolveMemory | Peak memory delta (bytes) |
| Status | OK/TIMEOUT/FAILED/SKIPPED |
| Residual | Infinity-norm of equation residual |

## Recommendations for future work

1. **Parallel DNFReduce branches**: `Or` branches in `DNFReduce` are independent and could be evaluated in parallel using `ParallelMap`. Note: parallel kernels do not share `$SolveCache`/`$ReduceCache`, so cache benefits would be per-branch only.
2. **Compiled Hamiltonian**: Use `Compile` or `FunctionCompile` for the Hamiltonian `H` and mass function `M` to reduce per-evaluation overhead (feasible for the default `alpha`, `g` but not user-overridden versions).
3. **Sparse matrix operations**: For large graphs (Grid1020), switch from dense `PseudoInverse` to sparse `LinearSolve` with pre-factorization.
