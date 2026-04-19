# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

**MFGraphs** is a Wolfram Language package for solving Mean Field Games on networks with congestion and switching costs. It converts network topology (vertices, edges, entry/exit flows, switching costs) into systems of equations, then solves for equilibrium flow distributions and value functions. Requires **Mathematica 12.0+**.

## Environment & Prerequisites

- **Mathematica 12.0 or later** — required for Wolfram Language syntax and `NDSolve`
- **wolframscript** — command-line tool included with Mathematica, used to run all scripts from the repository root
- **Paclet info**: See `PacletInfo.m` for package metadata (name: "MFGraphs", version: "0.0.2")

### Setting up wolframscript

Ensure `wolframscript` is on your PATH before running scripts:

```bash
# Check if available
which wolframscript

# If not found, add to PATH (edit ~/.zshrc, ~/.bash_profile, or ~/.bashrc):
export PATH="/Applications/Wolfram Desktop.app/Contents/MacOS:$PATH"
# or (for typical Mathematica installations):
export PATH="/opt/Wolfram/WolframEngine/12.0/Executables:$PATH"

# Verify installation and version
wolframscript -version
```

If `wolframscript` is not available, install Mathematica or Wolfram Engine from https://www.wolfram.com/download/.

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

### Built-in test cases

The package includes 34 pre-configured test networks accessible via `GetExampleData[key]`:

```mathematica
(* Access by numeric key *)
data = GetExampleData[12];                  (* Case 12 — 4-vertex attraction network *)
data = GetExampleData[7];                   (* Case 7 — solver stress test *)

(* Access by name (string) *)
data = GetExampleData["Paper example"];     (* Paper publication benchmark *)
data = GetExampleData["Grid0303"];          (* 3×3 grid network *)
data = GetExampleData["Jamaratv9"];         (* Large Jamarah network variant *)
data = GetExampleData["triangle with two exits"];  (* Symmetric stress case *)

(* Substitute symbolic parameters *)
data = GetExampleData[12] /. {I1 -> 100, U1 -> 0}  (* Replace symbolic parameters *)
```

For the complete list of 34 cases (both numeric and named), their definitions, and topologies, see `MFGraphs/Examples/ExamplesData.wl`. Many cases include symbolic parameters (e.g., `I1`, `U1`, `S1`) that must be substituted with numeric values before solving. Default substitutions for the benchmark suite are in `Scripts/BenchmarkHelpers.wls` (`$DefaultParams`).

### Testing a single case interactively

To verify package setup or debug a specific case, run in Mathematica/Wolfram Desktop:

```mathematica
<< MFGraphs`
data = GetExampleData[12] /. {I1 -> 100, U1 -> 0};
d2e = DataToEquations[data];
result = CriticalCongestionSolver[d2e];
IsFeasible[result]  (* should return True *)
```

### Quick Reference

Common operations at a glance:

| Task | Command |
|------|---------|
| Run fast tests (9 files, ~27 min) | `wolframscript -file Scripts/RunTests.wls fast` |
| Run all tests (19 files, slower) | `wolframscript -file Scripts/RunTests.wls all` |
| Run specific test file | `wolframscript -file MFGraphs/Tests/solver-contracts.mt` |
| Benchmark single case | `wolframscript -file Scripts/BenchmarkSuite.wls core case=12 timeout=300` |
| Regenerate API docs | `wolframscript -file Scripts/GenerateDocs.wls` |
| Clear solver cache | `ClearSolveCache[]` (in Mathematica) |
| Enable verbose output | `$MFGraphsVerbose = True` (in Mathematica) |

### Run tests and benchmarks

All scripts use `wolframscript` and must be run from the repository root.

#### Unit tests (correctness validation)

Fast regression suite for verifying package correctness:

```bash
# Fast suite (9 files, ~27 minutes)
wolframscript -file Scripts/RunTests.wls fast

# All tests (slower, more comprehensive)
wolframscript -file Scripts/RunTests.wls slow

# All test suites at once
wolframscript -file Scripts/RunTests.wls all
```

**Test execution expectations** (Mac with 8 cores, as of April 2026):
- `fast`: **~27 minutes** (9 test files: solver-contracts + monotone tests with full solver runs)
- `all`: **~60+ minutes** (all 19 test files, includes large-case solvers like Grid1010)
- Individual `.mt` file: **2–10 minutes** typical; `numeric-state.mt` and `critical-numeric-backend.mt` slower

**Test organization**:
- Test files are MUnit `.mt` files in `MFGraphs/Tests/` (19 total files)
- Run individual test files to isolate features:
  ```bash
  wolframscript -file MFGraphs/Tests/solver-contracts.mt       # Core solver validation
  wolframscript -file MFGraphs/Tests/numeric-state.mt          # Phase 5/6 Fictitious Play backend tests
  wolframscript -file MFGraphs/Tests/critical-numeric-backend.mt
  wolframscript -file MFGraphs/Tests/infeasible-status.mt      # Expected failure modes
  wolframscript -file MFGraphs/Tests/monotone-solver.mt        # Monotone operator solver
  ```
- Results go to standard Mathematica test output; no files are saved

#### Expected test behavior

Some test cases validate **failure modes** by design — these are not bugs:

- Tests in `infeasible-status.mt` verify that networks with insufficient capacity correctly return `"Infeasible"`
- Tests with switching costs that violate the triangle inequality are expected to fail feasibility validation (feature validation)
- Grid1020 (vlarge tier) is designed to hit `RecursionLimit` — intentional stress testing, not a regression

#### Performance benchmarks (tier-based)

Full benchmark suite across performance tiers:

```bash
# Full benchmark suite (all tiers)
wolframscript -file Scripts/BenchmarkSuite.wls

# Single tier (small/core/stress/paper/inconsistent-switching/large/vlarge)
wolframscript -file Scripts/BenchmarkSuite.wls small
wolframscript -file Scripts/BenchmarkSuite.wls core

# Single case with custom timeout
wolframscript -file Scripts/BenchmarkSuite.wls small case=1 timeout=60

# Track performance impact of changes (records in docs/history/DNF_PERFORMANCE_HISTORY.md)
wolframscript -file Scripts/BenchmarkSuite.wls --tag "description of change" --rationale "why" --changes "what changed"
```

**Benchmark tier specifications**:

| Tier | Cases | Timeout | Description |
|------|-------|---------|-------------|
| **small** | 1–6, 27 | 60s | Linear chains, simple cycles |
| **core** | 9, 10, 12, 14, 18 | 300s | Stable cases for routine benchmarking (recommended for CI) |
| **stress** | 7, 8, 104, "triangle with two exits" | 1800s | Symmetry-heavy or solver-stress cases |
| **paper** | "HRF Scenario 1" | 1800s | Paper publication benchmark reference |
| **inconsistent-switching** | 11, 15, 17, "Paper example" | 300s | Feature validation: switching costs violating triangle inequality |
| **large** | 13, 19–23, Braess variants, "Jamaratv9", Grid0303, Grid0404, Grid0505 | 900s | Multi-entrance/exit, larger graphs |
| **vlarge** | Grid0707, Grid0710, Grid1010, Grid1020 | 1800s | Stress test grids (up to 200 vertices, designed to hit RecursionLimit) |

**Benchmark notes**:
- Results (CSV, JSON, markdown) go to `Results/` (gitignored except `Results/compare_dnf.md`)
- Use `--tag` with `--rationale`, `--changes`, and `--interpretation` to track performance history in `docs/history/DNF_PERFORMANCE_HISTORY.md` and `docs/history/PARALLEL_PERFORMANCE_HISTORY.md`
- Test cases that fail validation due to switching costs violating the triangle inequality are **expected failures** (feature validation)
- `medium` is a legacy alias for `core` (supported for backward compatibility)

## Recent Changes & Highlights (April 2026)

**Test isolation and reliability improvements:**
- Fixed all remaining 9 test failures in critical congestion solver paths (commit 624136c)
- Implemented automatic test isolation to prevent state leakage between test cases
- Test runner now outputs failed test IDs directly for faster debugging
- Tests `J9F-critical_congestion.mt` and `Jm9-critical_congestion.mt` updated with proper resource cleanup

**Phase 5/6 Fictitious Play backend:**
- Numeric state machine fully implemented for cases where symbolic solver struggles
- 9 comprehensive tests in `numeric-state.mt` validate convergence and feasibility preservation
- **Not yet public API** — pending resolution of 4 technical debt issues (#72–#75)
- See Phase 5 & 6 section below for full details

## Architecture

### Pipeline overview

The package follows a linear pipeline: **network data → equations → solver**.

**Unified entrypoint** — `SolveMFG[data, Method -> m, opts]` accepts raw data or compiled equations and dispatches to the appropriate solver. With `Method -> "Automatic"`, it cascades: CriticalCongestion → Monotone → NonLinear, falling back on infeasibility or failure. When `"CongestionExponentFunction"` specifies alpha != 1 (scalar, Association, or Function), CriticalCongestion is skipped since it only handles alpha=1. Hamiltonian options (`"PotentialFunction"`, `"CongestionExponentFunction"`, `"InteractionFunction"`) and all downstream solver options are declared in `Options[SolveMFG]` and forwarded via `FilterRules`. `"CongestionExponentFunction"` accepts a scalar, `Function[edge, ...]`, or an `Association` mapping edges to exponents (unlisted edges default to 1); all forms are normalized by `NormalizeEdgeFunction`. Returns a standardized result envelope plus `"MethodUsed"` and `"MethodTrace"` keys.

**Individual solvers** can also be called directly:

```
Data (Association)
  ↓
DataToEquations[data]                    (DataToEquations.wl)
  ↓
CriticalCongestionSolver[d2e]           (Solvers.wl, uses DNFReduce.wl)
  ├─ DirectCriticalSolver[d2e]           (fast path for zero-switching-cost networks)
  └─ MFGSystemSolver + DNFSolveStep      (general case with switching costs)
  ↓
NonLinearSolver[criticalResult]  or  MonotoneSolverFromData[data]
```

**Solver selection**: `CriticalCongestionSolver` automatically dispatches to `DirectCriticalSolver` for networks with no switching costs (much faster); otherwise uses the full symbolic pipeline.

**MonotoneSolver vs MonotoneSolverFromData**: `MonotoneSolverFromData[data]` accepts raw data and calls `DataToEquations` internally. `MonotoneSolver[d2e]` accepts pre-compiled equations. `SolveMFG` uses `MonotoneSolver` when given compiled input, `MonotoneSolverFromData` otherwise.

### Module load order

Defined in `MFGraphs.wl`:
1. `Examples/ExamplesData.wl` — 34 built-in test cases via `GetExampleData[key]`
2. `Examples/TimeDependentExamples.wl` — Time-dependent example data
3. `DNFReduce.wl` — Boolean algebra solver (disjunctive normal form reduction with Solve/Reduce memoization, branch pruning, and post-reduction via `ReduceDisjuncts`/`SubsumptionPrune`)
4. `DataToEquations.wl` — Core converter: network topology → equations; implements `DataToEquations`, `MFGSystemSolver`, `MFGPreprocessing`, `TripleClean`
5. `Solvers.wl` — Critical-congestion solver suite extracted from `DataToEquations.wl` (Phase 3 refactor); implements `CriticalCongestionSolver`, `DirectCriticalSolver`, `IsCriticalSolution`
6. `NonLinearSolver.wl` — Iterative fixed-point solver (`NonLinearSolver`) using Hamiltonian framework (up to 15 iterations)
7. `Monotone.wl` — ODE-based gradient flow solver on Kirchhoff matrix using `NDSolve`
8. `TimeDependentSolver.wl` — Time-dependent MFG solver using backward-forward sweep on discretized time grid
9. `Graphics.wl` — Public visualization helpers: `NetworkGraphPlot`, `SolutionFlowPlot`, `ExitFlowPlot`

After submodule loading, `MFGraphs.wl` defines `SolveMFG` — the unified entrypoint (see Pipeline overview).

### Symbol contexts

Each submodule uses `Begin["`Private`"]` / `End[]` to isolate internal helpers:
- **Public symbols** (with `::usage` declarations): live in `MFGraphs`` context — these are the user-facing API
- **Internal symbols**: in `MFGraphs`Private`` — e.g., `$SolveCache`, `CachedSolve`, `TransitionsAt`, etc.

**Parallel dispatch gotcha**: When sending code to subkernels via `ParallelMap` or `ParallelTable`, `Function` is `HoldAll` (arguments aren't evaluated). Use `With[{val = expr}, Function[..., val ...]]` to inject values before dispatch.

### Inner solver pipeline

The critical congestion solver (in `Solvers.wl`, with supporting stages in `DataToEquations.wl`) uses a multi-stage pipeline:

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

All stationary solvers now return a **standardized result envelope** by default:

```mathematica
<| "Solver" -> "NonLinearSolver",
   "ResultKind" -> "Success"|"Failure"|"Degenerate"|"NonConverged",
   "Feasibility" -> "Feasible"|"Infeasible",
   "Message" -> "...",
   "Solution" -> {...},
   "ComparableFlowVector" -> {...},
   "KirchhoffResidual" -> ...,
   "AssoNonCritical" -> {...}  (* solver-specific payload key *)
|>
```

Solver-specific payload keys: `"AssoCritical"`, `"AssoNonCritical"`, `"AssoMonotone"`.

`SolveMFG` results additionally include `"MethodUsed"` and `"MethodTrace"`.

**Checking solutions**:
- `IsFeasible[result]` — accepts both legacy `"Status"` and standard `"Feasibility"`
- `IsNonLinearSolution[result]` — validates solution structure
- When `"Status"` is `"Infeasible"`, flow variables (`j[...]`) contain negative values (indicates no feasible solution)

### Choosing the right solver

| Solver | Best For | Speed | Notes |
|--------|----------|-------|-------|
| **CriticalCongestionSolver** | Symbolic, congestion exponent α=1, understanding structure | Seconds–minutes | Always try this first; serves as initial guess for other solvers. Fast path for zero-switching-cost networks. |
| **NonLinearSolver** | General non-linear (α≠1), fixed-point iteration | Minutes–hours | Requires potential function `V` defined. Starts from critical congestion seed. Up to 15 iterations. |
| **MonotoneSolver** | Gradient flow, ODE-based alternative | Minutes–hours | Uses Kirchhoff matrix. Good for comparison/validation. Slower but sometimes finds solutions NonLinear misses. |
| **SolveMFG with Method→"Automatic"** | Auto-routing, fallback cascades | Varies | Cascades: CriticalCongestion → Monotone → NonLinear. Adds `"MethodUsed"` and `"MethodTrace"` to results. |

**When you'll use Phase 5/6 (Fictitious Play backend):** If you add support for non-linear congestion exponents (α ≠ 1) to `CriticalCongestionSolver`, the symbolic `DNFReduce` pipeline may timeout or fail, triggering the numeric backend as fallback. Phase 5/6 prepares this handoff (currently internal-only pending issue resolution).

### Phase 5 & 6: Fictitious Play Backend (Internal v1)

**Status**: Implemented and tested, but **internal-only** (not yet wired into `CriticalCongestionSolver`). Kept v1 pending resolution of 4 technical debt issues (GitHub issues #72–#75).

**Phase 5: `SolveCriticalFictitiousPlayBackend`** (`MFGraphs.wl` line ~858)
- Iterative wrapper using `NestWhileList` that threads four Fictitious Play phases until `OracleReadyQ → True` or `MaxIterations` exhausted
- Pure functional state machine with immutable compound state Association
- Returns standardized backend result envelope with History, IterationLog, and full state snapshots
- Options: `MaxIterations` (20), `Temperature` (0.1), `Damping` (0.5), plus all Phase 4 thresholds
- **Use case**: Future support discovery for non-linear cases (α ≠ 1) when exact symbolic path fails

**Phase 6: `BuildOraclePrunedSystem`** (`DataToEquations.wl` line ~1922)
- Two-step Oracle pruning bridge translating `OracleState` from Phase 4 into a partially pruned `{EE, NN, OR}` triple
- Step 1: Reuses existing `BuildPrunedSystem` for OR-branch filtering
- Step 2: Injects explicit zero equalities for `PrunedZeroFlows` into EE
- Safe by design: only injects variables classified as inactive; strengthens constraints
- **Integration path (future)**: `OracleState` → `BuildOraclePrunedSystem` → pruned triple → `TripleClean` + `DNFSolveStep` → exact solver

**Testing**: 9 comprehensive tests in `MFGraphs/Tests/numeric-state.mt` (lines 539–681) covering convergence, feasibility preservation, history tracking, OR reduction, and end-to-end pipeline.

**Known limitations**:
- **DAG-only**: assumes acyclic topological order (violated for cyclic networks, causing instability)
- **Cost-scale fragility**: softmax transition probabilities poorly calibrated for heterogeneous switching costs
- **Float-to-exact chasm**: numeric flow output may lose precision on >5-digit costs when converting to symbolic
- **Cycle recovery**: no mechanism to detect and recover support when DAG assumption breaks

**Future work**: Address vulnerabilities (#72–#75) before public API exposure. Wire Phase 5 as fallback in `CriticalCongestionSolver` when α ≠ 1 forces numeric-to-symbolic bridge.

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

Performance history tracked in `docs/history/DNF_PERFORMANCE_HISTORY.md`. When modifying `DNFReduce.wl`, run `CompareDNF.wls --tag "description"` to record impact.

### Parallelization

Parallelized operations use `ParallelMap`/`ParallelTable` when workload exceeds `$MFGraphsParallelThreshold` (default 6). Kernels launch lazily (only when a parallel path is actually taken).

Parallelized sites: `Simplify` over switching-cost inequalities, `Select`/`Simplify` over transition flows, `Reduce` over Or-branches in `DNFSolveStep`, and `SolveMassRoot` grid in `PrecomputeM`.

**Parallel gotchas**:
- `ParallelNeeds["pkg`"]` is a no-op if context exists; use `DistributeDefinitions` instead
- `Block[{sym}]` clears all `DownValues` — provide explicit defaults when wrapping solver parameters

Performance history tracked in `docs/history/PARALLEL_PERFORMANCE_HISTORY.md`. Use `BenchmarkSuite.wls --tag "description"` to record entries.

## Utility Scripts

The `Scripts/` directory contains specialized tools for benchmarking, profiling, and performance analysis. All scripts use `wolframscript` and must be run from the repository root.

### Performance comparison and tracking

**`CompareDNF.wls`** — Compare DNF solver performance before/after code changes:
```bash
# Establish baseline before making changes
wolframscript -file Scripts/CompareDNF.wls --tag "baseline"

# Make code changes...

# Measure impact after changes
wolframscript -file Scripts/CompareDNF.wls --tag "after optimization"

# Compare performance on a single case
wolframscript -file Scripts/CompareDNF.wls case=7
```
Results are recorded in `docs/history/DNF_PERFORMANCE_HISTORY.md` and `Results/compare_dnf.md`.

**`CompareSolvers.wls`** — Compare performance across all three solvers (CriticalCongestion, NonLinear, Monotone) on representative test cases.

### Profiling and bottleneck analysis

**`BottleneckReport.wls`** — Generate detailed timing breakdown per solver function:
```bash
wolframscript -file Scripts/BottleneckReport.wls
```
Produces `Results/bottleneck_report.md` with call counts and timing for each profiled stage (preprocessing, DNFReduce, TripleClean, etc.). Use this to identify performance bottlenecks.

**`ProfilePreprocessing.wls`** — Deep-dive profiling of the preprocessing stage (`MFGPreprocessing` + `TripleClean`).

**`ProfileReductionStrategies.wls`** — Compare different disjunct reduction strategies (`ReduceDisjuncts` vs `SubsumptionPrune`) on performance-critical cases.

**`ProfileInstrument.wls`** — Generate instrumentation data for custom analysis.

### Helper scripts

**`BenchmarkHelpers.wls`** — Loaded by `BenchmarkSuite.wls`; defines tier configurations, default parameters, and helper functions. Contains `$DefaultParams` substitutions and tier definitions.

**`RunBenchmarkCI.wls`** — Runs `small` and `core` tiers in CI mode (in-process isolation), writes results to `Results/benchmark_latest.json`. Use this for quick automated validation across the stable case set.

**`GenerateDocs.wls`** — Regenerate auto-generated documentation (e.g., `API_REFERENCE.md` from `::usage` declarations).

**`GenerateReferenceHashes.wls`** — Generate reference solution hashes for verification.

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

Results are recorded in `docs/history/DNF_PERFORMANCE_HISTORY.md` and `Results/compare_dnf.md`.

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

## Code style conventions

- Avoid single-letter `System` symbols (e.g., `K`, `D`, `C`) as variable names — they shadow built-ins. Use descriptive names or suffixes like `KM`.
- For CLI scripts in `Scripts/`, suppress Wolfram Code Analysis warnings for `Print[]` with inline lint blocks:
  ```wolfram
  (* :!CodeAnalysis::BeginBlock:: *)
  (* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)
  ...
  (* :!CodeAnalysis::EndBlock:: *)
  ```
- `API_REFERENCE.md` is auto-generated — edit `::usage` strings in source, then run `wolframscript -file Scripts/GenerateDocs.wls`
- `docs/history/DNF_PERFORMANCE_HISTORY.md` and `docs/history/PARALLEL_PERFORMANCE_HISTORY.md` are auto-appended by benchmark scripts — do not manually edit

## CI pipeline

None. GitHub Actions workflows were removed in commit `64a3ee9` — there is no `.github/workflows/` directory. All test and benchmark runs must be performed locally per the "Run tests and benchmarks" section above before opening a PR.

## Git workflow

**Never commit directly to `master`.** All changes must be made on a new branch and merged via pull request.
