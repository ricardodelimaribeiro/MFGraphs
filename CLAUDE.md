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

The two-phase load is intentional: `Scenario.wl` and `ExampleScenarios.wl` must be available before `BeginPackage` because submodules depend on their exports.

Archived/inactive modules include:
- `MFGraphs/Examples/archive/ExamplesData.wl`
- `MFGraphs/archive/DataToEquations.wl`

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

Run tests from repository root:

```bash
wolframscript -file Scripts/RunTests.wls fast
wolframscript -file MFGraphs/Tests/scenario-kernel.mt
wolframscript -file MFGraphs/Tests/make-unknowns.mt
```

Run a single `.mt` file (note: `Scripts/RunSingleTest.wls` has a hardcoded `$RepoRoot` — update it if your repo is not at `/Users/ribeirrd/Documents/GitHub/MFGraphs`):

```bash
wolframscript -file Scripts/RunSingleTest.wls MFGraphs/Tests/scenario-kernel.mt
```

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
