# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

**MFGraphs** is a Wolfram Language package for solving Mean Field Games on networks with congestion and switching costs. It converts network topology (vertices, edges, entry/exit flows, switching costs) into systems of equations, then solves for equilibrium flow distributions and value functions.

## Running tests and benchmarks

All scripts use `wolframscript`:

```bash
# Run full benchmark suite (validates all 31+ test cases)
wolframscript -file Scripts/BenchmarkSuite.wls

# Run only a specific tier: small, medium, large, or vlarge
wolframscript -file Scripts/BenchmarkSuite.wls small

# Profile solver bottlenecks
wolframscript -file Scripts/BottleneckReport.wls

# Compare DNF performance across git commits, with optional tag
wolframscript -file Scripts/CompareDNF.wls --tag "Description of change"
```

Test tiers by case complexity:
- **small**: Cases 1–6, 27 (linear chains, cycle) — 60s timeout
- **medium**: Cases 7–15, 17–18, 104, triangle — 300s timeout
- **large**: Cases 13, 19–23, Braess variants, Jamaratv9, Grid0303 — 900s timeout
- **vlarge**: Grid1020 — 1800s timeout (hits RecursionLimit by design)

Test cases that fail validation due to switching costs not satisfying the triangle inequality are expected failures (the feature is working correctly).

## Architecture

The package has a linear pipeline: **network data → equations → solver**.

```
Data (Association)
  ↓
Data2Equations[data]          (D2E2.m)
  ↓
CriticalCongestionSolver      (D2E2.m, uses DNFReduce.m)
  ↓
NonLinear / MonotoneSolverFromData
```

**Module load order** (defined in `MFGraphs.m`):
1. `Examples/ExamplesData.m` — 34 built-in test cases via `DataG[key]`
2. `DNFReduce.m` — Boolean algebra solver (disjunctive normal form reduction with Solve/Reduce memoization and branch pruning)
3. `D2E2.m` — Core converter: network topology → equations; implements `Data2Equations`, `CriticalCongestionSolver`, `MFGSystemSolver`, `TripleClean`
4. `NonLinearSolver.m` — Iterative fixed-point solver using Hamiltonian framework (up to 15 iterations)
5. `Monotone.m` — ODE-based gradient flow solver on Kirchhoff matrix using `NDSolve`

**Legacy files** (`D2E2-currents.m`, `IterationFunction2.m`) are kept for reference but **not loaded** by the package.

## Solver outputs

- `CriticalCongestionSolver` returns an association with key `"AssoCritical"` — the zero-flow equilibrium
- `NonLinear` returns an association with key `"AssoNonCritical"` — the general congestion solution
- Check solution validity with `IsNonLinearSolution[result]`

## Key solver parameters (set before calling `NonLinear`)

```mathematica
$MFGraphsVerbose = False     (* suppress timing/progress messages *)
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

Symbolic parameters (e.g., `I1`, `U1`, `S1`) can be used and substituted with `/.` rules. The `DataG[key]` function returns networks with symbolic parameters that must be substituted before solving.

## DNFReduce performance

`DNFReduce` is the primary solver bottleneck. It uses:
- **Solver-aware memoization**: caches `Solve`/`Reduce` calls to avoid redundant symbolic computation
- **Branch pruning**: prunes false branches early to avoid exponential blowup

Performance history is tracked in `DNF_PERFORMANCE_HISTORY.md`. When making changes to `DNFReduce.m`, run `CompareDNF.wls` with `--tag` to record the impact.
