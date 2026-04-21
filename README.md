# MFGraphs

**MFGraphs** is a Wolfram Language package for studying **Mean Field Games on networks** with congestion and switching costs. The package uses an undirected core graph representation with directed-flow variables, converts network data into a standardized symbolic system, solves the active runtime problem class, and provides plotting helpers for inspecting the resulting flows.

At the **`v0.4.0-critical-only`** milestone, the repository is intentionally narrowed to the **critical-congestion solver surface**. The active package path now focuses on the critical specialization with edgewise `alpha = 1`, while legacy non-critical solver families are preserved only as archive tags for historical retrieval.

## Quick navigation

- [Repository layout](#repository-layout)
- [Installation](#installation)
- [Status at `v0.4.0-critical-only`](#status-at-v040-critical-only)
- [Quick start](#quick-start)
- [Defining a network](#defining-a-network)
- [Solvers and result objects](#solvers-and-result-objects)
- [Plotting](#plotting)
- [Scenario API](#scenario-api)
- [Built-in examples](#built-in-examples)
- [Running tests](#running-tests)
- [Package structure](#package-structure)
- [Archive tags for legacy solvers](#archive-tags-for-legacy-solvers)
- [Documentation](#documentation)
- [License](#license)

For contributor workflow notes, see [CLAUDE.md](CLAUDE.md). For troubleshooting, see [TROUBLESHOOTING.md](TROUBLESHOOTING.md).

## Repository layout

The repository is organized around a small runtime package, focused scripts, and supporting research and planning material.

| Path | Purpose |
|---|---|
| `MFGraphs/` | Package source, examples, tests, and kernel initialization |
| `Scripts/` | Test runners, benchmarks, comparisons, and documentation utilities |
| `Results/` | Benchmark and profiling outputs, mostly generated artifacts |
| `docs/` | Internal planning notes, validation logs, and history files |
| `repro/` | Targeted reproduction and verification scripts |
| `research/` | Reference papers and related supporting material |

## Installation

Clone the repository and load the package from Mathematica or WolframScript.

```mathematica
(* Option 1: if the package directory is already on $Path *)
Needs["MFGraphs`"]

(* Option 2: direct load from a checkout *)
Get["/path/to/MFGraphs/MFGraphs/MFGraphs.wl"]
```

The package is intended for Mathematica 12.0 or later.

## Status at `v0.4.0-critical-only`

This tag marks the point where the runtime package surface was deliberately simplified and synchronized around the **critical-congestion workflow**. In practical terms, `MFGraphs\`` now loads the critical solver path, the public graphics helpers, and the scenario module, while non-critical solver families are no longer part of the active package load path.

| Area | Status at `v0.4.0-critical-only` |
|---|---|
| Solver scope | **Critical-only** active runtime surface |
| Main solver entrypoints | `SolveMFG`, `CriticalCongestionSolver` |
| Data compilation | `DataToEquations` |
| Public plotting API | `NetworkGraphPlot`, `SolutionFlowPlot`, `ExitFlowPlot` |
| Scenario support | `makeScenario`, `validateScenario`, `completeScenario`, `scenarioQ`, `ScenarioData` |
| Legacy non-critical solvers | Removed from active runtime path and preserved through archive tags |

## Quick start

The simplest workflow is to load example data, compile it with `DataToEquations`, and solve it with the critical solver.

```mathematica
<< MFGraphs`

Data = GetExampleData[12] /. {I1 -> 100, U1 -> 0};
d2e = DataToEquations[Data];
result = CriticalCongestionSolver[d2e];

result["Feasibility"]
result["Solution"]
```

If you want a single entrypoint, `SolveMFG` accepts either raw model data or a compiled `DataToEquations` association.

```mathematica
result1 = SolveMFG[Data];
result2 = SolveMFG[d2e];
```

Solver outputs are standardized associations. In typical use, the most important keys are the feasibility classification, a message, and the recovered solution data.

| Common key | Meaning |
|---|---|
| `"Solver"` | Solver/backend name used for the run |
| `"ResultKind"` | High-level result classification |
| `"Feasibility"` | Feasibility verdict such as `"Feasible"` or `"Infeasible"` |
| `"Message"` | Short diagnostic summary |
| `"Solution"` | Primary solution association |
| `"ComparableFlowVector"` | Numeric flow representation for comparisons |
| `"KirchhoffResidual"` | Residual measure for flow conservation checks |

## Defining a network

Raw MFGraphs model data is represented as an `Association`. The package expects a network topology, entry flows, exit costs, and switching-cost data.

```mathematica
Data = <|
  "Vertices List" -> {1, 2, 3, 4},
  "Adjacency Matrix" -> {
    {0, 1, 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, 1},
    {0, 0, 0, 0}
  },
  "Entrance Vertices and Flows" -> {{1, 200}},
  "Exit Vertices and Terminal Costs" -> {{4, 0}},
  "Switching Costs" -> {{1, 2, 3, 5}, {3, 2, 1, 5}}
|>;
```

| Key | Description |
|---|---|
| `"Vertices List"` | Vertex labels used throughout the model |
| `"Adjacency Matrix"` | Square 0-1 matrix encoding the directed edges |
| `"Entrance Vertices and Flows"` | `{{vertex, flow}, ...}` pairs for inflow conditions |
| `"Exit Vertices and Terminal Costs"` | `{{vertex, cost}, ...}` pairs for terminal conditions |
| `"Switching Costs"` | `{{from, at, to, cost}, ...}` edge-switch penalties; use `{}` when absent |

Symbolic parameters such as `I1`, `U1`, or `S1` may be left symbolic and substituted later with replacement rules.

## Solvers and result objects

At this tag, the active solver story is intentionally narrow. The package is designed around the **critical-congestion** case, and the exported solver entrypoints are documented with that scope in mind.

```mathematica
Data = GetExampleData[7] /. {I1 -> 50, U1 -> 0, U2 -> 0};
d2e = DataToEquations[Data];
result = CriticalCongestionSolver[d2e];

IsFeasible[result]
result["AssoCritical"]
```

`CriticalCongestionSolver` operates on the standardized `DataToEquations` output. `SolveMFG` provides a convenience wrapper that routes through the same critical-only runtime surface. Utility routines such as `MFGPreprocessing`, `MFGSystemSolver`, `IsCriticalSolution`, and flow-feasibility helpers remain available for lower-level analysis and debugging.

## Plotting

The first wave of the public graphics API is extracted into `MFGraphs/Graphics.wl`. These helpers are intended to cover the common tasks of visualizing the network topology, solved edge flows, and exit-flow totals.

```mathematica
result = CriticalCongestionSolver[d2e];

NetworkGraphPlot[d2e]
SolutionFlowPlot[d2e, result["Solution"]]
```

| Function | Purpose |
|---|---|
| `NetworkGraphPlot[d2e]` | Draws the network structure from the compiled model |
| `SolutionFlowPlot[d2e, solution]` | Draws the network with edge styling and labels derived from solved net flows |
| `ExitFlowPlot[exitFlows]` | Produces a bar chart for total flow exiting at each exit vertex |

Additional workbook-oriented graphics helpers, such as density and value-function plots, are still candidates for later extraction.

## Scenario API

The package now also includes a lightweight typed scenario layer. This is useful when you want to package a raw model, parameter substitutions, metadata, and validation state into a single structured object.

```mathematica
sc = makeScenario[<|
  "Identity" -> <|"name" -> "demo"|>,
  "Model" -> Data,
  "Data" -> {I1 -> 100, U1 -> 0}
|>];

scenarioQ[sc]
ScenarioData[sc, "Identity"]
```

| Function | Purpose |
|---|---|
| `makeScenario[assoc]` | Validates, completes, and wraps raw scenario data |
| `validateScenario[scenario]` | Checks required top-level and model-level structure |
| `completeScenario[scenario]` | Fills derived metadata such as content hash and benchmark defaults |
| `scenarioQ[x]` | Predicate for typed scenario objects |
| `ScenarioData[scenario, key]` | Accessor for scenario contents |

## Built-in examples

Use `GetExampleData[key]` to retrieve predefined benchmark and teaching examples. These are stored in `MFGraphs/Examples/ExamplesData.wl` and remain the quickest way to exercise the package.

| Key | Network |
|---|---|
| `7`, `8` | Y-network with one entrance and two exits, without or with switching |
| `11`, `12` | Attraction problem, without or with switching |
| `"Braess congest"` | Braess-paradox congestion example |
| `"Jamaratv9"` | Jamarat pilgrimage network |
| `"Paper example"` | Four-vertex example with two entrances, two exits, and switching |

## Running tests

For routine regression checks, use the centralized test runner.

```bash
wolframscript -file Scripts/RunTests.wls fast
wolframscript -file Scripts/RunTests.wls slow
```

The current test surface includes package-loading checks, graphics API regression coverage, and critical-only solver regression cases, including the inconsistent-switching critical recovery behavior added during the `v0.4.0-critical-only` cleanup.

## Package structure

The runtime package is now relatively compact. The loader pulls in examples, symbolic reduction, data compilation, the critical solver path, plotting helpers, and the scenario module.

```text
MFGraphs/
  MFGraphs.wl             Main package loader and exported public API
  DNFReduce.wl            DNF reduction and symbolic logical simplification
  DataToEquations.wl      Raw network data -> standardized equation association
  Solvers.wl              Critical-congestion solver suite
  Graphics.wl             Public plotting API
  Scenario.wl             Typed scenario kernel
  Examples/
    ExamplesData.wl       Built-in network examples
  Tests/
    *.mt                  MUnit regression tests
  Kernel/
    init.m                Paclet initialization
Scripts/
  RunTests.wls            Test-suite runner
  BenchmarkSuite.wls      Benchmark driver
  CompareDNF.wls          DNF comparison utility
  GenerateDocs.wls        API documentation generator
```

## Archive tags for legacy solvers

Legacy non-critical solver files are **not** part of the active runtime path at `v0.4.0-critical-only`, but they remain recoverable through archive tags in the Git history.

| Ref | Purpose |
|---|---|
| `v0.4.0-critical-only` | Current milestone tag for the critical-only package surface |
| `archive/non-critical-solvers` | Archived non-critical solver state |
| `archive/non-critical-solvers-full` | Fuller archived non-critical solver state |
| `archive/nonlinear-residual-bug` | Historical diagnostic tag related to nonlinear residual behavior |

Useful retrieval commands are shown below.

```bash
# Check out the critical-only milestone
git checkout v0.4.0-critical-only

# Inspect or recover archived legacy files
git show archive/non-critical-solvers:MFGraphs/NonLinearSolver.wl
git checkout archive/non-critical-solvers -- MFGraphs/NonLinearSolver.wl
git show archive/non-critical-solvers-full:MFGraphs/Monotone.wl
git checkout archive/non-critical-solvers-full -- MFGraphs/Monotone.wl
```

## Documentation

The generated API reference is available in [API_REFERENCE.md](API_REFERENCE.md). That document is derived from the package `::usage` strings and is the best reference for the currently exported public symbols.

## License

This project is part of ongoing research. Please contact the authors before using it in publications.
