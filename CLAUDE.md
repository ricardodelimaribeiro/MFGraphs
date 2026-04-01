# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

**MFGraphs** is a Wolfram Language package for solving Mean Field Games on networks with congestion and switching costs. It converts network topology (vertices, edges, entry/exit flows, switching costs) into systems of equations, then solves for equilibrium flow distributions and value functions. Requires **Mathematica 12.0+**.

## Environment & Prerequisites

- **Mathematica 12.0 or later** — required for Wolfram Language syntax and `NDSolve`
- **wolframscript** — command-line tool included with Mathematica, used to run all scripts from the repository root
- **Optional**: GitHub Actions setup requires `WOLFRAM_CI_ID` and `WOLFRAM_CI_PASS` environment secrets for running tests via `.github/workflows/tests.yml`
- **Paclet info**: See `PacletInfo.m` for package metadata (name: "MFGraphs", version: "0.0.2")

## Quick Start

### Load the package

```mathematica
(* In Mathematica *)
Get["path/to/MFGraphs/MFGraphs.wl"]
(* or use Needs if paclet is registered: *)
Needs["MFGraphs`"]
```

### Solve a simple example

```mathematica
(* Get a built-in test case *)
data = GetExampleData["Y_switch"];
d2e = DataToEquations[data];
result = CriticalCongestionSolver[d2e];

(* Check if feasible and get solution *)
IsFeasible[result]
If[IsFeasible[result], Print["Flows: ", result["AssoCritical"]]]
```

### Run tests and benchmarks

All scripts use `wolframscript` and must be run from the repository root:

```bash
# Fast regression suite (9 files, ~27 minutes)
wolframscript -file Scripts/RunTests.wls fast

# All tests (slower, more comprehensive)
wolframscript -file Scripts/RunTests.wls slow

# Full benchmark suite across all tiers
wolframscript -file Scripts/BenchmarkSuite.wls

# Single tier (small/medium/large/vlarge)
wolframscript -file Scripts/BenchmarkSuite.wls small

# Track performance impact of changes
wolframscript -file Scripts/CompareDNF.wls --tag "description of change"
```

**Test tier specifications**:

| Tier | Cases | Timeout | Notes |
|------|-------|---------|-------|
| **small** | 1–6, 27 | 60s | Linear chains, cycles |
| **medium** | 7–15, 17–18, 104, triangle, Paper | 300s | Mid-size networks |
| **large** | 13, 19–23, Braess variants, Jamaratv9, Grid0303 | 900s | Large networks |
| **vlarge** | Grid1020 | 1800s | Stress test (RecursionLimit by design) |

Expected times on a typical Mac (as of March 2025):
- `RunTests.wls fast`: **~27 minutes** (solver-contracts + monotone tests with full solver runs)
- `RunTests.wls slow`: **significantly longer** (includes large-case solvers)

**Test infrastructure notes**:
- Benchmark results (CSV, JSON, markdown) go to `Results/` (gitignored except `Results/compare_dnf.md`)
- Test cases that fail validation due to switching costs violating the triangle inequality are **expected failures** (feature validation)
- Run benchmarks with `--tag` to track performance history; see `DNF_PERFORMANCE_HISTORY.md` and `PARALLEL_PERFORMANCE_HISTORY.md`

## Architecture

### Pipeline overview

The package follows a linear pipeline: **network data → equations → solver**.

```
Data (Association)
  ↓
DataToEquations[data]                    (DataToEquations.wl)
  ↓
CriticalCongestionSolver[d2e]           (DataToEquations.wl, uses DNFReduce.wl)
  ↓
NonLinearSolver[criticalResult]  or  MonotoneSolverFromData[data]
```

### Module load order

Defined in `MFGraphs.wl`:
1. `Examples/ExamplesData.wl` — 34 built-in test cases via `GetExampleData[key]`
2. `DNFReduce.wl` — Boolean algebra solver (disjunctive normal form reduction with Solve/Reduce memoization, branch pruning, and post-reduction via `ReduceDisjuncts`/`SubsumptionPrune`)
3. `DataToEquations.wl` — Core converter: network topology → equations; implements `DataToEquations`, `CriticalCongestionSolver`, `MFGSystemSolver`, `TripleClean`
4. `NonLinearSolver.wl` — Iterative fixed-point solver (`NonLinearSolver`) using Hamiltonian framework (up to 15 iterations)
5. `Monotone.wl` — ODE-based gradient flow solver on Kirchhoff matrix using `NDSolve`
6. `TimeDependentSolver.wl` — Time-dependent MFG solver using backward-forward sweep on discretized time grid

### Symbol contexts

Each submodule uses `Begin["`Private`"]` / `End[]` to isolate internal helpers:
- **Public symbols** (with `::usage` declarations): live in `MFGraphs`` context — these are the user-facing API
- **Internal symbols**: in `MFGraphs`Private`` — e.g., `$SolveCache`, `CachedSolve`, `TransitionsAt`, etc.

**Parallel dispatch gotcha**: When sending code to subkernels via `ParallelMap` or `ParallelTable`, `Function` is `HoldAll` (arguments aren't evaluated). Use `With[{val = expr}, Function[..., val ...]]` to inject values before dispatch.

### Inner solver pipeline

The critical congestion solver inside `DataToEquations.wl` uses a multi-stage pipeline:

```
CriticalCongestionSolver
  → MFGPreprocessing     — builds "InitRules" (partial solution) and "NewSystem" {EE, NN, OR}
    → TripleClean        — fixed-point loop of TripleStep: solve equalities, substitute, repeat
  → MFGSystemSolver      — applies flow approximation, then:
    → TripleClean        — simplify again with substituted flows
    → DNFSolveStep       — calls DNFReduce on the Or-branch, then TripleClean on result
```

System decomposition (triple `{EE, NN, OR}`):
- **EE** — equalities (solved and substituted by `TripleStep`)
- **NN** — inequalities (non-negative flows, switching bounds)
- **OR** — disjunctions (complementarity conditions from optimality)

`TripleClean` is `FixedPoint[TripleStep, ...]` — repeatedly solves equalities and substitutes until no new equalities emerge.

### Solver chain & return formats

`NonLinearSolver` expects the output of `CriticalCongestionSolver` (reads `"AssoCritical"` key). Benchmark suite follows: `DataToEquations → CriticalCongestionSolver → NonLinearSolver`.

**Legacy format** (default):
- `CriticalCongestionSolver`: returns `"AssoCritical"` (zero-flow equilibrium) + `"Status"` (`"Feasible"` or `"Infeasible"`)
- `NonLinearSolver`: returns `"AssoNonCritical"` (general congestion solution) + `"Status"`
- `MonotoneSolverFromData`: returns bare solution, `Null`, or message association

**Standard format** (add `"ReturnShape" -> "Standard"`):
```mathematica
result = NonLinearSolver[d2e, "ReturnShape" -> "Standard"]
(* Returns normalized envelope: *)
<| "Solver" -> "NonLinearSolver",
   "ResultKind" -> "...",
   "Feasibility" -> "Feasible"|"Infeasible",
   "Message" -> "...",
   "Solution" -> {...},
   "AssoNonCritical" -> {...}  (* solver-specific payload *)
|>
```

**Checking solutions**:
- `IsFeasible[result]` — accepts both legacy `"Status"` and standard `"Feasibility"`
- `IsNonLinearSolution[result]` — validates solution structure
- When `"Status"` is `"Infeasible"`, flow variables (`j[...]`) contain negative values (indicates no feasible solution)

## Configuration

Set these global parameters before calling solvers to customize behavior:

```mathematica
$MFGraphsVerbose = False              (* enable/disable timing and progress messages *)
$MFGraphsParallelThreshold = 6        (* minimum list length to trigger parallel dispatch *)
$MFGraphsParallelReady = False        (* set True after DistributeDefinitions to subkernels *)
$ReduceDisjunctsThreshold = 200       (* disjunct count cutoff: Full Reduce vs Subsumption *)

alpha[edge_] := 1                     (* congestion exponent *)
g[m_, edge_] := -1/m^2               (* interaction potential *)
V = Function[{x, edge}, 0.5 Sin[2 Pi (x + 1/4)]^2]  (* potential function *)
```

**Cache management**:
- `ClearSolveCache[]` — clears memoization between different problem instances (important to prevent stale results when solving multiple unrelated networks)

## Network data format

```mathematica
Data = <|
  "Vertices List"                    -> {1, 2, 3, 4},
  "Adjacency Matrix"                 -> {{0,1,0,0}, {0,0,1,0}, {0,0,0,1}, {0,0,0,0}},
  "Entrance Vertices and Flows"      -> {{1, 200}},
  "Exit Vertices and Terminal Costs" -> {{4, 0}},
  "Switching Costs"                  -> {{1,2,3, 5}, {3,2,1, 5}}
|>;
```

Symbolic parameters (e.g., `I1`, `U1`, `S1`) can be used and substituted with `/.` rules. `GetExampleData[key]` returns networks with symbolic parameters that must be substituted before solving. Default substitutions for benchmarks are in `Scripts/BenchmarkHelpers.wls` (`$DefaultParams`).

## Performance optimization notes

### DNFReduce bottleneck

`DNFReduce` is the primary solver bottleneck. Optimizations include:
- **Solver-aware memoization**: `CachedSolve`/`CachedReduce` hash inputs (SHA256) to avoid redundant symbolic computation
- **Branch pruning**: prunes false branches early to avoid exponential blowup
- **And-Or early exit**: stops immediately if any Or-branch leaves `xp` unchanged (means `xp` already satisfies the disjunction)
- **Post-reduction via `ReduceDisjuncts`**: removes redundant branches. For small output (≤ `$ReduceDisjunctsThreshold`, default 200) uses `Reduce[Or @@ branches, Reals]`; for larger output applies `SubsumptionPrune` (pairwise implication checks) first, then `Reduce`

**Important**: Call `ClearSolveCache[]` between different problem instances to prevent stale cached results.

Performance history tracked in `DNF_PERFORMANCE_HISTORY.md`. When modifying `DNFReduce.wl`, run `CompareDNF.wls --tag "description"` to record impact.

### Parallelization

Parallelized operations use `ParallelMap`/`ParallelTable` when workload exceeds `$MFGraphsParallelThreshold` (default 6). Kernels launch lazily (only when a parallel path is actually taken).

Parallelized sites: `Simplify` over switching-cost inequalities, `Select`/`Simplify` over transition flows, `Reduce` over Or-branches in `DNFSolveStep`, and `SolveMassRoot` grid in `PrecomputeM`.

**Parallel gotchas**:
- `ParallelNeeds["pkg`"]` is a no-op if context exists; use `DistributeDefinitions` instead
- `Block[{sym}]` clears all `DownValues` — provide explicit defaults when wrapping solver parameters

Performance history tracked in `PARALLEL_PERFORMANCE_HISTORY.md`. Use `BenchmarkSuite.wls --tag "description"` to record entries.

## Debugging & Profiling Workflows

### Enable verbose logging

```mathematica
$MFGraphsVerbose = True;
result = CriticalCongestionSolver[d2e];  (* prints timing and progress *)
```

### Debug infeasible solutions

```mathematica
result = CriticalCongestionSolver[d2e];
If[!IsFeasible[result],
  Print["Infeasibility detected. Negative flows:"];
  Print[Select[result["AssoCritical"], # < 0 &]]
]
```

### Profile solver bottlenecks

```bash
# Run profiling on representative cases from each tier
wolframscript -file Scripts/BottleneckReport.wls
```

This generates detailed timing breakdowns for each solver stage (preprocessing, DNFReduce, TripleClean, etc.).

### Track performance regressions

```bash
# Before making changes to DNFReduce.wl or solvers, establish baseline
wolframscript -file Scripts/CompareDNF.wls --tag "before optimization"

# Make changes, then measure impact
wolframscript -file Scripts/CompareDNF.wls --tag "after optimization"

# Compare a single test case
wolframscript -file Scripts/CompareDNF.wls Y_switch
```

Results are recorded in `DNF_PERFORMANCE_HISTORY.md` and `Results/compare_dnf.md`.

### Debug parallel kernel issues

```mathematica
(* Ensure parallel kernels are initialized *)
EnsureParallelKernels[];

(* Verify definitions are distributed *)
DistributeDefinitions[alpha, g, V];
$MFGraphsParallelReady = True;

(* Run solver and check parallel execution *)
$MFGraphsVerbose = True;
result = NonLinearSolver[d2e];
```

If parallel dispatch doesn't trigger, check:
- List length exceeds `$MFGraphsParallelThreshold` (default 6)
- Relevant definitions are sent to subkernels via `DistributeDefinitions`
- `$MFGraphsParallelReady = True` is set

### Clear solver cache between problem instances

```mathematica
(* Solve first problem *)
result1 = CriticalCongestionSolver[d2e1];

(* Must clear cache before solving a different problem *)
ClearSolveCache[];

(* Solve second problem *)
result2 = CriticalCongestionSolver[d2e2];
```

Failure to clear the cache can cause stale memoized results to leak into unrelated problems.

## Git workflow

**Never commit directly to `master`.** All changes must be made on a new branch and merged via pull request.
