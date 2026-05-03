# MFGraphs

**MFGraphs** is a Wolfram Language package for Mean Field Games on networks.

The repository is currently in a **core scenario-kernel phase**. The active package surface is centered on:
- typed scenario construction
- example scenario factories
- exact symbolic unknown bundle construction
- structural system construction
- critical-congestion symbolic solving
- visualization helpers

## Current status

Loaded by `Needs["MFGraphs`"]`:
- `MFGraphs/primitives.wl`
- `MFGraphs/scenarioTools.wl`
- `MFGraphs/examples.wl`
- `MFGraphs/unknownsTools.wl`
- `MFGraphs/systemTools.wl`
- `MFGraphs/solversTools.wl`
- `MFGraphs/orchestrationTools.wl`
- `MFGraphs/graphicsTools.wl`

Archived/inactive modules live under `MFGraphs/archive/`.

## Solver Status

The symbolic solvers (`MFGraphs/solversTools.wl`) are designed for **critical congestion only** (`Alpha = 1` on every edge). Non-critical systems (`Alpha != 1` or edge-specific non-1 `EdgeAlpha`) fail explicitly. The default user-facing path is `solveScenario`, which uses `dnfReduceSystem`; raw `reduceSystem` remains available as a direct Wolfram `Reduce` baseline.

Hamiltonian parameters `V`, `G`, `EdgeV`, and `EdgeG` are validated and preserved on scenarios for future density-per-edge visualization work, but current structural system construction applies only `Alpha` and `EdgeAlpha`.

## Installation

```mathematica
Needs["MFGraphs`"]
(* or *)
Get["/path/to/MFGraphs/MFGraphs/MFGraphs.wl"]
```

Requires Mathematica / Wolfram Language 12.0+.

## Quick start (core phase)

```mathematica
Needs["MFGraphs`"];

s = makeScenario[<|
  "Model" -> <|
    "Vertices" -> {1, 2, 3},
    "Adjacency" -> {{0, 1, 0}, {0, 0, 1}, {0, 0, 0}},
    "Entries" -> {{1, 10}},
    "Exits" -> {{3, 0}},
    "Switching" -> {}
  |>
|>];

unk = makeSymbolicUnknowns[s];
sys = makeSystem[s, unk];
sol = solveScenario[s];

scenarioQ[s]
symbolicUnknownsQ[unk]
mfgSystemQ[sys]
isValidSystemSolution[sys, sol]
```

Example factory workflow:

```mathematica
s = getExampleScenario[12, {{1, 100}}, {{4, 0}}];
```

## Active API (current phase)

Scenario kernel:
- `scenario`, `scenarioQ`
- `makeScenario`, `validateScenario`, `completeScenario`
- `scenarioData`

Example scenario factories:
- `getExampleScenario`
- `gridScenario`, `cycleScenario`, `graphScenario`, `amScenario`

Symbolic unknown/system kernels:
- `symbolicUnknowns`, `symbolicUnknownsQ`, `symbolicUnknownsData`, `makeSymbolicUnknowns`
- `mfgSystem`, `mfgSystemQ`, `systemData`, `makeSystem`
- `buildBoundaryData`, `buildFlowData`, `buildComplementarityData`, `buildHamiltonianData`
- `solveScenario`, `dnfReduceSystem`, `reduceSystem`, `isValidSystemSolution`

Core symbolic primitives:
- `j`, `u`, `z`, `alpha`, `Cost`

## Intentionally unavailable in this phase

Not loaded by current `MFGraphs.wl` bootstrap:
- `ScenarioByKey`, `GetExampleData`
- `DataToEquations`, `CriticalCongestionSolver`
- legacy/extended solver modules and solver benchmarking paths

## Running tests

From repository root:

```bash
wolframscript -file Scripts/RunTests.wls fast
wolframscript -file MFGraphs/Tests/scenario-kernel.mt
wolframscript -file MFGraphs/Tests/symbolic-unknowns.mt
```

Current active runner suites (`Scripts/RunTests.wls`):
- `fast`: `scenario-kernel.mt`, `symbolic-unknowns.mt`, `reduce-system.mt`, `scenario-consistency.mt`, `graphicsTools.mt`, `orchestration.mt`
- `all`: alias for `fast`
- `archive`: archived compatibility/legacy suites (explicit use only)
- `full`: `fast + archive`

## Repository structure

```text
MFGraphs/
  MFGraphs.wl
  primitives.wl
  scenarioTools.wl
  examples.wl
  unknownsTools.wl
  systemTools.wl
  solversTools.wl
  orchestrationTools.wl
  graphicsTools.wl
  archive/
  Tests/
    scenario-kernel.mt
    symbolic-unknowns.mt
    archive/
  Kernel/init.m
Scripts/
  RunTests.wls
```

## Documentation

- Canonical workflow/reference: [CLAUDE.md](CLAUDE.md)
- Pointer companion: [GEMINI.md](GEMINI.md)
- Assistant pointer policy: [AGENTS.md](AGENTS.md)
- Troubleshooting: [TROUBLESHOOTING.md](TROUBLESHOOTING.md)

## License

This project is part of ongoing research. Please contact the authors before using it in publications.
