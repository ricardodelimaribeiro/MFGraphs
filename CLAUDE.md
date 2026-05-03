# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

**MFGraphs** is currently in a **core scenario-kernel phase**.
The active package surface focuses on:
- typed scenario construction (`scenarioTools.wl`)
- example scenario factories (`examples.wl`)
- exact symbolic unknown bundle construction (`unknownsTools.wl`)
- structural system construction (`systemTools.wl`)
- symbolic solver (`solversTools.wl`)
- shared primitives (`primitives.wl`)

## Environment & prerequisites

- **Mathematica / Wolfram Language 12.0+**
- **wolframscript** available on PATH
- **Local Wolfram app/docs path**: `/Applications/Wolfram.app`
- Package metadata in `PacletInfo.m` (`Name -> "MFGraphs"`, `Version -> "0.0.2"`)

Quick check:

```bash
which wolframscript
wolframscript -version
```

## Current package load behavior

`Needs["MFGraphs`"]` loads the full core stack via the numerics-style architecture:

```wolfram
PrependTo[$Path, DirectoryName[$InputFileName]];
BeginPackage["MFGraphs`", {"primitives`", "scenarioTools`", "examples`",
                            "unknownsTools`", "systemTools`", "solversTools`",
                            "orchestrationTools`", "graphicsTools`"}];
EndPackage[]
```

`EndPackage` appends all dependency contexts to `$ContextPath`, so every public symbol is accessible unqualified after a single `Needs["MFGraphs`"]`.

Loading DAG (each file is an independent package with its own flat context):
```
primitives`  →  scenarioTools`  →  examples`
                               →  unknownsTools`  →  systemTools`  →  solversTools`
                                                                   →  orchestrationTools`
                                                                   →  graphicsTools`
```

Archived/inactive modules live under `MFGraphs/archive/`.

## Core architecture

### Composition pipeline

Every workflow follows the same three-step chain:

```
makeScenario[input]  →  makeSymbolicUnknowns[s]  →  makeSystem[s, unk]
```

- **`makeScenario`** normalizes input (Graph objects → adjacency matrix), validates, fills defaults, then builds and caches `buildAuxiliaryTopology` output inside the returned `scenario[<|...|>]` wrapper.
- **`makeSymbolicUnknowns`** reads the cached topology and generates symbolic variable families (`j`, `u`, `z`).
- **`makeSystem`** orchestrates four typed builders — `buildBoundaryData`, `buildFlowData`, `buildComplementarityData`, `buildHamiltonianData` — each returning its own typed record (`mfgBoundaryData`, etc.), merged into the final `mfgSystem`.

### Typed wrapper pattern

All kernel objects follow one pattern: `TypeHead[Association]` with a predicate and accessor:

| Head | Predicate | Accessor |
|---|---|---|
| `scenario` | `scenarioQ` | `scenarioData[s]` / `scenarioData[s, key]` |
| `symbolicUnknowns` | `symbolicUnknownsQ` | `symbolicUnknownsData[u]` / `symbolicUnknownsData[u, key]` |
| `mfgSystem` | `mfgSystemQ` | `systemData[sys]` / `systemData[sys, key]` |

`typeData[obj, key]` returns `Missing["KeyAbsent", key]` for absent keys. System sub-records (`mfgBoundaryData` etc.) are typed records nested inside `mfgSystem`; use `systemData[sys, key]` or `systemDataFlatten[sys]` for flattened access to their fields.

### Topology construction and caching

`buildAuxiliaryTopology` creates synthetic auxiliary vertices as string labels (`"auxEntry1"`, `"auxEntry2"`, …) for network boundary points and computes all valid three-vertex transitions (triples `{vIn, vMid, vOut}`) for flow variable generation. Building topology is expensive — it is computed once inside `makeScenario` and cached in `scenarioData[s, "Topology"]`. Always pass the full scenario object to downstream constructors rather than rebuilding.

Edges are network connections, not directed movement permissions. A non-symmetric AM is symmetrized (`Unitize[AM + Transpose[AM]]`) before topology construction, so `makeScenario` with AM or `AM + Transpose[AM]` produces identical topologies. Movement restrictions are imposed via switching costs = `Infinity` on the blocked triples.

### Switching costs

Switching costs may be supplied as a List of 4-tuples or an Association with 3-tuple keys. Values may be numeric or `Infinity` for blocked transitions. `completeScenario` normalizes all representations and fills every missing triple with cost `0`; if called without cached topology, it warns and rebuilds topology from the model before completing.

### Examples registry

`getExampleScenario[key]` returns a 6-argument `Function[{entries, exits, sc, alpha, V, g}, ...]`. Integer keys 1–23 plus named strings ("Braess split", "Jamaratv9", etc.) are registered in `$ExampleScenarios`. Call with 3 arguments to use canonical defaults: `getExampleScenario[key, entries, exits]`.

## Active public API (current phase)

### Scenario kernel (`scenarioTools``)
- `scenario`, `scenarioQ`
- `makeScenario`, `validateScenario`, `completeScenario`
- `scenarioData`
- `buildAuxiliaryTopology`, `deriveAuxPairs`, `buildAuxTriples`

### Example scenario factories (`examples``)
- `getExampleScenario`
- `gridScenario`, `cycleScenario`, `graphScenario`, `amScenario`

### Symbolic unknown/system kernels (`unknownsTools``, `systemTools``)
- `symbolicUnknowns`, `symbolicUnknownsQ`, `symbolicUnknownsData`, `makeSymbolicUnknowns`
- `mfgSystem`, `mfgSystemQ`, `systemData`, `systemDataFlatten`, `makeSystem`
- **System Builders**: `buildBoundaryData`, `buildFlowData`, `buildComplementarityData`, `buildHamiltonianData`
- **System Records**: `mfgBoundaryData`, `mfgFlowData`, `mfgComplementarityData`, `mfgHamiltonianData`
- **Linear Helpers**: `getKirchhoffLinearSystem`, `getKirchhoffMatrix`

### Solver (`solversTools`)
All public solver entrypoints share a common preprocessing pipeline (`buildSolverInputs`) and only support critical congestion systems (`Alpha == 1` on every edge).

- Default user-facing orchestration uses `solveScenario`, which calls `dnfReduceSystem`.
- `reduceSystem` — calls `Reduce[constraints, allVars, Reals]` directly
- `dnfReduceSystem` — equality-substitution + disjunction-distribution via `dnfReduce` (avoids `Reduce` timeouts)
- `optimizedDNFReduceSystem` — opt-in exact branch-state DNF reducer for small branch families, with exact DNF fallback for larger residual variable sets
- `activeSetReduceSystem` — opt-in exact active-set branch reducer for the critical-congestion complementarity structure, with exact DNF fallback for larger residual variable sets
- `booleanReduceSystem` — converts to DNF via `BooleanConvert`, then calls `Reduce` independently per disjunct; options: `"DisjunctTimeout"` (default 30s), `"ReturnAll"` (default `False`)
- `findInstanceSystem` — calls `FindInstance[constraints, allVars, Reals]` after accumulated linear preprocessing; option: `"Timeout"` (default `Infinity`)
- `dnfReduce` — internal simplifier used by `dnfReduceSystem`; eliminates equalities and distributes over disjunctions
- `isValidSystemSolution` — solution validator; option `"ReturnReport" -> True` gives per-block breakdown

### Solver design notes

- **Common preprocessing pipeline** (`buildSolverInputs`): all solvers call this first. It collects structural equations, applies exit-value rules, and runs `accumulateLinearRules` to eliminate linear equalities symbolically. `activeSetReduceSystem` bypasses `buildSolverInputs` and instead splits constraints manually into deterministic (balance, boundary, inequalities) and active (complementarity) groups, accumulating linear rules on the deterministic block only before processing the complementarity structure.

- **Branch-state variable threshold**: `optimizedDNFReduceSystem` and `activeSetReduceSystem` choose their branch-state path only when `Length[allVars] <= 8`. Above this threshold both fall back to the standard `dnfReduce` path. This threshold is a fixed internal constant, not an option.

- **`branchStateReduce` vs `dnfReduce`**: `dnfReduce` is a recursive, expression-level simplifier that eliminates equalities and distributes over disjunctions on the fly, building up a DNF accumulator. `branchStateReduce` is a fold-based branch state machine — it maintains an explicit list of `<|"Rules" -> ..., "Residuals" -> ...|>` associations, one per surviving branch, and adds each conjunct to every live branch, pruning branches that become `False`. The two approaches can produce different residual shapes and branch orderings.

- **`dnfReduceDiagnosticReport` global state**: the instrumented reducer writes into `$dnfCurrentMonitor` via mutable `=` assignments during evaluation. This is intentionally non-reentrant — only one diagnostic report should run at a time. Interrupted calls leave the global dirty; this is acceptable for script/test use but would be unsafe in a concurrent or nested context.

- **`solutionResultKind` taxonomy**: classifies solver output into `"Rules"` (fully determined list of rules, or association with `Equations -> True`), `"NoSolution"` (empty list or `Equations -> False`), `"Branched"` (residual containing `Or`), `"Parametric"` (tracked `j`/`u` variables remain free but no disjunction), `"Residual"` (other symbolic residual), `"TIMEOUT"`, or `"Unknown"`. Used by benchmark scripts and profilers for result reporting.

- **Topology sensitivity**: grid cases with the same vertex/edge count can differ significantly in solver time depending on layout orientation. The 2×3 and 3×2 grids both have 6 vertices and 7 edges, but `dnfReduceSystem` takes ~4.4× longer on grid-2x3 than grid-3x2 (see BENCHMARKS.md). This reflects the number of complementarity alternatives the DNF path must enumerate.

### Orchestration (`orchestrationTools`)
- `solveScenario` — builds symbolic unknowns and system, then solves with `dnfReduceSystem` by default
- `SolveMFG` — compatibility wrapper for raw association input, delegated through `solveScenario`

### Graphics (`graphicsTools`)
- `scenarioTopologyPlot` — entry/exit/internal vertex coloring
- `mfgSolutionPlot` — coordinated grid of flow, value, and density views
- `mfgFlowPlot` — flow-only network plot with directed real and auxiliary edges
- `mfgValuePlot` — real-edge endpoint and interpolated interior value labels
- `mfgDensityPlot` — real-edge inferred agent density labels and styling
- `mfgTransitionPlot` — transition graph (nodes = edges, edges = j[a,b,c])
- `augmentAuxiliaryGraph` — road-traffic graph derived from AuxPairs/AuxTriples
- `mfgAugmentedPlot` — paper infrastructure graph (nodes = (e,v) pairs)

### Runtime flags

- `$MFGraphsVerbose` — set `True` to enable progress/timing prints (default `False`)
- `$MFGraphsParallelThreshold` — minimum list length for `ParallelMap`/`ParallelTable` (default `6`; set `Infinity` to disable)

### Core symbolic primitives
- `j`, `u`, `z`, `alpha`, `Cost`

## Intentionally unavailable in this phase

The following are currently not loaded from `MFGraphs.wl`:
- `ScenarioByKey`, `GetExampleData`
- `DataToEquations`, `CriticalCongestionSolver`
- legacy/extended solver modules and solver benchmarking workflows
- solver-phase backend helpers (archived in `MFGraphs/archive/SolverBackendHelpers.wl`)

**V, G, EdgeV, EdgeG Hamiltonian parameters** are validated and stored in the scenario
schema but are **not applied** in `buildHamiltonianData` yet — only `Alpha`/`EdgeAlpha` are
used in current system construction. They are reserved for future density-per-edge
visualization work. Current `EqGeneral` is the integrated transport-equation form
associated with the Hamiltonian equation: the integral over `[0,1]` of the adjoint
of the linearization of the stationary Hamilton-Jacobi equation.

If solver work resumes, restore the module loading and docs in a dedicated pass.

## Quick start (core phase)

```mathematica
Needs["MFGraphs`"];

s = makeScenario[<|
  "Model" -> <|
    "Vertices"  -> {1, 2, 3},
    "Adjacency" -> {{0,1,0},{0,0,1},{0,0,0}},
    "Entries"   -> {{1, 10}},
    "Exits"     -> {{3, 0}},
    "Switching" -> {}
  |>
|>];

unk = makeSymbolicUnknowns[s];
sys = makeSystem[s, unk];
sol = solveScenario[s];
optSol = solveScenario[s, optimizedDNFReduceSystem];
activeSol = solveScenario[s, activeSetReduceSystem];

scenarioQ[s]
symbolicUnknownsQ[unk]
mfgSystemQ[sys]
isValidSystemSolution[sys, sol]
```

Example factory usage:

```mathematica
s = getExampleScenario[12, {{1, 100}}, {{4, 0}}];
```

## Testing

Active suite (`Scripts/RunTests.wls`):
- `fast`: `scenario-kernel.mt`, `symbolic-unknowns.mt`, `reduce-system.mt`, `scenario-consistency.mt`, `graphicsTools.mt`, `orchestration.mt`
- `all`: alias for `fast`
- `archive`: archived compatibility/legacy suites (explicit use only)
- `full`: `fast + archive`

Run tests from repository root using the current runner workflow:

```bash
# Run the active fast suite
wolframscript -file Scripts/RunTests.wls fast

# Run a single test file
wolframscript -file Scripts/RunSingleTest.wls MFGraphs/Tests/reduce-system.mt
```

`all` currently resolves to the active core-phase set:

```bash
wolframscript -file Scripts/RunTests.wls all
```

Solver-oriented suites are intentionally excluded in this phase while solver
modules remain out of scope.

## Benchmarking

```bash
# Benchmark the DNF-first default solver across representative cases (results → Results/)
wolframscript -file Scripts/BenchmarkSystemSolver.wls

# With options: tag appends an entry to BENCHMARKS.md; timeout sets per-case limit
wolframscript -file Scripts/BenchmarkSystemSolver.wls --tag "after my change" --timeout 60

# Compare opt-in exact solvers
wolframscript -file Scripts/BenchmarkSystemSolver.wls --solver optimizeddnf
wolframscript -file Scripts/BenchmarkSystemSolver.wls --solver activeset

# Run a single bench case
wolframscript -file Scripts/BenchmarkSystemSolver.wls --case chain-3v-1exit

# Non-mutating staged profiler (construction + solver timing per stage)
wolframscript -file Scripts/ProfileScenarioKernel.wls --case example-12 --timeout 10
wolframscript -file Scripts/ProfileScenarioKernel.wls --case example-12 --solver optimizeddnf --timeout 10

# DNF reducer diagnostic profiler (branch/ordering analysis, does not write results)
wolframscript -file Scripts/ProfileDNFReduce.wls --case grid-3x2 --order original
wolframscript -file Scripts/ProfileDNFReduce.wls --case all --order all --timeout 60

# Raw Reduce baseline benchmark, retained for comparison
wolframscript -file Scripts/BenchmarkReduceSystem.wls
```

Available bench cases: `chain-2v`, `chain-3v-1exit`, `chain-3v-2exit`, `example-7`, `chain-5v-1exit`, `grid-2x3`, `grid-3x2`, `example-12`.

Benchmark results land in `Results/` as timestamped CSV + WL solution files and as latest files. Solutions store the full symbolic output per case alongside timing. `ProfileScenarioKernel.wls` prints diagnostics to stdout and does not write results.

## Docs generation

`API_REFERENCE.md` is auto-generated from `::usage` strings — do not edit it directly:

```bash
wolframscript -file Scripts/GenerateDocs.wls
```

## Development notebook

`MFGraphs/Getting started.wl` is an interactive notebook script — not intended to be run end-to-end with `wolframscript`. It includes force-reload logic for safe iteration and demonstrates the full Scenario → Unknowns → System pipeline with formatted output. Use it as a sandbox; run it cell-by-cell in the Mathematica front end.

## Documentation maintenance policy

- Treat this file as the **canonical** project workflow/reference doc.
- Keep `GEMINI.md` as a short pointer doc to avoid duplicated guidance.
- When workflow/API changes, update this file first.

## Code style notes

- Avoid single-letter variable names that collide with WL/System symbols.
- Prefer typed constructors/accessors (`makeScenario`, `scenarioData`, `makeSymbolicUnknowns`, `makeSystem`) over ad-hoc association wiring.
- Keep tests deterministic and context-safe; prefer qualified symbol checks (e.g. `scenarioTools`scenarioData`) for package-surface assertions.

## Git workflow

- Never commit directly to `master`.
- Create feature branches and merge through pull requests.
