# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

**MFGraphs** is a Wolfram Language package for solving Mean Field Games on networks with congestion and switching costs. It converts network topology (vertices, edges, entry/exit flows, switching costs) into systems of equations, then solves for equilibrium flow distributions and value functions. Requires Mathematica 12.0+.

## Running tests and benchmarks

All scripts use `wolframscript` and must be run from the repository root:

```bash
# Run full benchmark suite (validates all 31+ test cases)
wolframscript -file Scripts/BenchmarkSuite.wls

# Run curated regression suites
wolframscript -file Scripts/RunTests.wls fast
wolframscript -file Scripts/RunTests.wls slow

# Run only a specific tier: small, medium, large, or vlarge
wolframscript -file Scripts/BenchmarkSuite.wls small

# Profile solver bottlenecks (runs representative cases from each tier)
wolframscript -file Scripts/BottleneckReport.wls

# Compare DNF performance across git commits, with optional tag
wolframscript -file Scripts/CompareDNF.wls --tag "Description of change"

# Compare a single case
wolframscript -file Scripts/CompareDNF.wls Y_switch
```

Test tiers by case complexity:
- **small**: Cases 1–6, 27 (linear chains, cycle) — 60s timeout
- **medium**: Cases 7–15, 17–18, 104, triangle, Paper — 300s timeout
- **large**: Cases 13, 19–23, Braess variants, Jamaratv9, Grid0303 — 900s timeout
- **vlarge**: Grid1020 — 1800s timeout (hits RecursionLimit by design)

Expected wall-clock times (on a typical Mac, as of March 2025):
- `RunTests.wls fast` (38 tests across 9 files): **~27 minutes** — includes solver-contracts and monotone tests which run full CriticalCongestionSolver, NonLinearSolver, and MonotoneSolver (NDSolve) calls
- `RunTests.wls slow` (3 files): significantly longer due to large-case solvers

Test cases that fail validation due to switching costs not satisfying the triangle inequality are expected failures (the feature is working correctly).

Benchmark results go to `Results/` (CSV, JSON, markdown). These generated files are gitignored except for `Results/compare_dnf.md`.

## Architecture

The package has a linear pipeline: **network data → equations → solver**.

```
Data (Association)
  ↓
DataToEquations[data]                    (DataToEquations.wl)
  ↓
CriticalCongestionSolver[d2e]           (DataToEquations.wl, uses DNFReduce.wl)
  ↓
NonLinearSolver[criticalResult]  or  MonotoneSolverFromData[data]
```

**Module load order** (defined in `MFGraphs.wl`):
1. `Examples/ExamplesData.wl` — 34 built-in test cases via `GetExampleData[key]`
2. `DNFReduce.wl` — Boolean algebra solver (disjunctive normal form reduction with Solve/Reduce memoization and branch pruning)
3. `DataToEquations.wl` — Core converter: network topology → equations; implements `DataToEquations`, `CriticalCongestionSolver`, `MFGSystemSolver`, `TripleClean`
4. `NonLinearSolver.wl` — Iterative fixed-point solver (`NonLinearSolver`) using Hamiltonian framework (up to 15 iterations)
5. `Monotone.wl` — ODE-based gradient flow solver on Kirchhoff matrix using `NDSolve`

Each submodule uses `Begin["`Private`"]` / `End[]` to isolate internal helpers. Public symbols (those with `::usage`) are declared before `Begin` and live in the `MFGraphs`` context. Internal symbols like `$SolveCache`, `CachedSolve`, `TransitionsAt`, etc. are in `MFGraphs`Private``.

### Inner solver pipeline (DataToEquations.wl)

The critical congestion solver has a multi-stage pipeline inside `DataToEquations.wl`:

```
CriticalCongestionSolver
  → MFGPreprocessing     — builds "InitRules" (partial solution) and "NewSystem" {EE, NN, OR}
    → TripleClean        — fixed-point loop of TripleStep: solve equalities, substitute, repeat
  → MFGSystemSolver      — applies flow approximation, then:
    → TripleClean        — simplify again with substituted flows
    → DNFSolveStep          — calls DNFReduce on the Or-branch, then TripleClean on result
```

The system is decomposed into a triple `{EE, NN, OR}` by `SystemToTriple`:
- **EE** — equalities (solved and substituted by `TripleStep`)
- **NN** — inequalities (non-negative flows, switching bounds)
- **OR** — disjunctions (complementarity conditions from optimality)

`TripleClean` is `FixedPoint[TripleStep, ...]` — it repeatedly solves equalities and substitutes until no new equalities emerge.

### Solver chain

`NonLinearSolver` expects the output of `CriticalCongestionSolver` as input (it reads the `"AssoCritical"` key). The benchmark suite follows this chain: `DataToEquations → CriticalCongestionSolver → NonLinearSolver`.

## Solver outputs

- By default, solvers keep their legacy return format.
- `CriticalCongestionSolver` legacy output includes `"AssoCritical"` (zero-flow equilibrium) and `"Status"` (`"Feasible"` or `"Infeasible"`)
- `NonLinearSolver` legacy output includes `"AssoNonCritical"` (general congestion solution) and `"Status"`
- `MonotoneSolver` / `MonotoneSolverFromData` legacy output is a bare solution association, `Null`, or a degenerate message association
- Add `"ReturnShape" -> "Standard"` to any solver call to receive a normalized envelope with `"Solver"`, `"ResultKind"`, `"Feasibility"`, `"Message"`, and `"Solution"`
- Standard results retain the solver-specific payload key when available: `"AssoCritical"`, `"AssoNonCritical"`, or `"AssoMonotone"`
- Check feasibility with `IsFeasible[result]` — it accepts both legacy `"Status"` and standardized `"Feasibility"`
- Check solution validity with `IsNonLinearSolution[result]`
- When `"Status"` is `"Infeasible"`, flow variables (`j[...]`) contain negative values indicating the problem has no feasible non-negative flow solution

## Key solver parameters (set before calling `NonLinearSolver`)

```mathematica
$MFGraphsVerbose = False     (* suppress timing/progress messages — default is False *)
alpha[edge_] := 1            (* congestion exponent *)
g[m_, edge_] := -1/m^2      (* interaction potential *)
V = Function[{x, edge}, 0.5 Sin[2 Pi (x + 1/4)]^2]  (* potential function — user-defined *)
```

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

Symbolic parameters (e.g., `I1`, `U1`, `S1`) can be used and substituted with `/.` rules. The `GetExampleData[key]` function returns networks with symbolic parameters that must be substituted before solving. Default substitutions for benchmarks are in `Scripts/BenchmarkHelpers.wls` (`$DefaultParams`).

## DNFReduce performance

`DNFReduce` is the primary solver bottleneck. It uses:
- **Solver-aware memoization**: `CachedSolve`/`CachedReduce` hash inputs (SHA256) to avoid redundant symbolic computation
- **Branch pruning**: prunes false branches early to avoid exponential blowup
- **And-Or early exit**: when distributing over Or-branches sequentially, stops immediately if any branch leaves `xp` unchanged (meaning `xp` already satisfies the disjunction)

**Important**: Call `ClearSolveCache[]` between different problem instances to prevent stale cached results.

Performance history is tracked in `DNF_PERFORMANCE_HISTORY.md`. When making changes to `DNFReduce.wl`, run `CompareDNF.wls` with `--tag` to record the impact.

## Git workflow

**Never commit directly to `master`.** All changes must be made on a new branch and merged via pull request.
