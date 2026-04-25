# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

**MFGraphs** is currently in a **core scenario-kernel phase**.
The active package surface focuses on:
- typed scenario construction (`Scenario.wl`)
- example scenario factories (`Examples/ExampleScenarios.wl`)
- unknown bundle construction (`Unknowns.wl`)
- structural system construction (`System.wl`)

Solver modules are intentionally not loaded in this phase.

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

`Needs["MFGraphs`"]` loads the current core stack:
- Preloaded via `Get` before `BeginPackage`: `MFGraphs/Scenario.wl`, `MFGraphs/Examples/ExampleScenarios.wl` (each establishes `MFGraphs`` independently)
- Loaded via `Scan[Get, ...]` inside `Private` after `BeginPackage`: `MFGraphs/Unknowns.wl`, `MFGraphs/System.wl`

The two-phase load is intentional: `ExampleScenarios.wl` calls `makeScenario` during factory construction, so `Scenario.wl` must be fully evaluated before `MFGraphs.wl`'s own `BeginPackage` runs. Moving either Phase 1 module into Phase 2 (or vice versa) breaks cross-module symbol resolution.

Archived/inactive modules include:
- `MFGraphs/Examples/archive/ExamplesData.wl`
- `MFGraphs/archive/DataToEquations.wl`

## Core architecture

### Composition pipeline

Every workflow follows the same three-step chain:

```
makeScenario[input]  →  makeUnknowns[s]  →  makeSystem[s, unk]
```

- **`makeScenario`** normalizes input (Graph objects → adjacency matrix), validates, fills defaults, then builds and caches `BuildAuxiliaryTopology` output inside the returned `scenario[<|...|>]` wrapper.
- **`makeUnknowns`** reads the cached topology and generates symbolic variable families (`j`, `u`, `z`).
- **`makeSystem`** orchestrates four typed builders — `BuildBoundaryData`, `BuildFlowData`, `BuildComplementarityData`, `BuildHamiltonianData` — each returning its own typed record (`mfgBoundaryData`, etc.), merged into the final `mfgSystem`.

### Typed wrapper pattern

All kernel objects follow one pattern: `TypeHead[Association]` with a predicate and accessor:

| Head | Predicate | Accessor |
|---|---|---|
| `scenario` | `scenarioQ` | `ScenarioData[s]` / `ScenarioData[s, key]` |
| `unknowns` | `unknownsQ` | `UnknownsData[u]` / `UnknownsData[u, key]` |
| `mfgSystem` | `mfgSystemQ` | `SystemData[sys]` / `SystemData[sys, key]` |

`TypeData[obj, key]` returns `Missing["KeyAbsent", key]` for absent keys. Sub-records (`mfgBoundaryData` etc.) follow the same accessor convention.

### Topology construction and caching

`BuildAuxiliaryTopology` creates synthetic auxiliary vertices as string labels (`"auxEntry1"`, `"auxEntry2"`, …) for network boundary points and computes all valid three-vertex transitions (triples `{vIn, vMid, vOut}`) for flow variable generation. Building topology is expensive — it is computed once inside `makeScenario` and cached in `ScenarioData[s, "Topology"]`. Always pass the full scenario object to downstream constructors rather than rebuilding.

The graph is always undirected. A non-symmetric AM is symmetrized (`Unitize[AM + Transpose[AM]]`) before topology construction, so `makeScenario` with AM or `AM + Transpose[AM]` produces identical topologies. Edge directionality is imposed exclusively via switching costs = `Infinity` on the blocked triples.

### Switching costs

Switching costs may be supplied as a List of 4-tuples or an Association with 3-tuple keys. `completeScenario` normalizes all representations and fills every missing triple with cost `0`. A triangle-inequality check runs at validation time; violations emit a message but do **not** abort scenario construction.

### Examples registry

`GetExampleScenario[key]` returns a 6-argument `Function[{entries, exits, sc, alpha, V, g}, ...]`. Integer keys 1–23 plus named strings ("Braess split", "Jamaratv9", etc.) are registered in `$ExampleScenarios`. Call with 3 arguments to use canonical defaults: `GetExampleScenario[key, entries, exits]`.

## Active public API (current phase)

### Scenario kernel
- `scenario`, `scenarioQ`
- `makeScenario`, `validateScenario`, `completeScenario`
- `ScenarioData`

### Example scenario factories
- `GetExampleScenario`
- `GridScenario`, `CycleScenario`, `GraphScenario`, `AMScenario`

### Unknowns/system kernels
- `unknowns`, `unknownsQ`, `UnknownsData`, `makeUnknowns`
- `mfgSystem`, `mfgSystemQ`, `SystemData`, `SystemDataFlatten`, `makeSystem`
- **System Builders**: `BuildBoundaryData`, `BuildFlowData`, `BuildComplementarityData`, `BuildHamiltonianData`
- **System Records**: `mfgBoundaryData`, `mfgFlowData`, `mfgComplementarityData`, `mfgHamiltonianData`
- **Linear Helpers**: `GetKirchhoffLinearSystem`, `GetKirchhoffMatrix`

### Shared topology helpers (from `Scenario.wl`)
- `BuildAuxiliaryTopology`, `DeriveAuxPairs`, `BuildAuxTriples`

### Runtime flags
- `$MFGraphsVerbose` — set `True` to enable progress/timing prints (default `False`)
- `$MFGraphsParallelThreshold` — minimum list length for `ParallelMap`/`ParallelTable` (default `6`; set `Infinity` to disable)

### Core symbolic primitives
- `j`, `u`, `z`, `alpha`, `Cost`

## Intentionally unavailable in this phase

The following are currently not loaded from `MFGraphs.wl`:
- `ScenarioByKey`, `GetExampleData`
- `DataToEquations`, `CriticalCongestionSolver`, `SolveMFG`
- legacy/extended solver modules and solver benchmarking workflows
- solver-phase backend helpers (archived in `MFGraphs/archive/SolverBackendHelpers.wl`)

**V, G, EdgeV, EdgeG Hamiltonian parameters** are validated and stored in the scenario
schema but are **not applied** in `BuildHamiltonianData` — only `Alpha`/`EdgeAlpha` are
used in system construction. Implementing V and G requires extending the HJ-equation
builder. Treat these fields as schema-only for now.

If solver work resumes, restore the module loading and docs in a dedicated pass.

## Quick start (core phase)

```mathematica
Needs["MFGraphs`"];

s = makeScenario[<|
  "Model" -> <|
    "Vertices List" -> {1, 2, 3},
    "Adjacency Matrix" -> {{0,1,0},{0,0,1},{0,0,0}},
    "Entrance Vertices and Flows" -> {{1, 10}},
    "Exit Vertices and Terminal Costs" -> {{3, 0}},
    "Switching Costs" -> {}
  |>
|>];

unk = makeUnknowns[s];
sys = makeSystem[s, unk];

scenarioQ[s]
unknownsQ[unk]
mfgSystemQ[sys]
```

Example factory usage:

```mathematica
f = GetExampleScenario[12];
s = f[{{1, 100}}, {{4, 0}}, {}, 1, 0, Function[z, -1/z]];
```

## Testing

Active suite (`Scripts/RunTests.wls`):
- `fast`: `scenario-kernel.mt`, `make-unknowns.mt`
- `slow`, `legacy-fast`, `legacy-slow`: empty in current phase

Run tests from repository root using the current runner workflow:

```bash
wolframscript -file Scripts/RunTests.wls fast
```

`all` currently resolves to the same core-phase set (no active slow suites):

```bash
wolframscript -file Scripts/RunTests.wls all
```

Solver-oriented suites are intentionally excluded in this phase while solver
modules remain out of scope.

## Development notebook

`MFGraphs/Getting started.wl` is an interactive notebook script — not intended to be run end-to-end with `wolframscript`. It includes force-reload logic for safe iteration and demonstrates the full Scenario → Unknowns → System pipeline with formatted output. Use it as a sandbox; run it cell-by-cell in the Mathematica front end.

## Documentation maintenance policy

- Treat this file as the **canonical** project workflow/reference doc.
- Keep `GEMINI.md` as a short pointer doc to avoid duplicated guidance.
- When workflow/API changes, update this file first.

## Code style notes

- Avoid single-letter variable names that collide with WL/System symbols.
- Prefer typed constructors/accessors (`makeScenario`, `ScenarioData`, `makeUnknowns`, `makeSystem`) over ad-hoc association wiring.
- Keep tests deterministic and context-safe; prefer qualified symbol checks (`MFGraphs` context) for package-surface assertions.

## Git workflow

- Never commit directly to `master`.
- Create feature branches and merge through pull requests.
