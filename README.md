# MFGraphs

**MFGraphs** is a Wolfram Language package for Mean Field Games on networks.

The repository is currently in a **core scenario-kernel phase**. The active package surface is centered on:
- typed scenario construction
- example scenario factories
- unknown-variable bundle construction
- structural system construction

Solver modules are intentionally not loaded by the current package bootstrap.

## Current status

Loaded by `Needs["MFGraphs`"]`:
- `MFGraphs/Scenario.wl`
- `MFGraphs/Examples/ExampleScenarios.wl`
- `MFGraphs/Unknowns.wl`
- `MFGraphs/System.wl`

Archived/inactive modules:
- `MFGraphs/Examples/archive/ExamplesData.wl`
- `MFGraphs/archive/DataToEquations.wl`

## Solver Status

The symbolic solver (`MFGraphs/solver.wl`) is currently designed for **critical congestion only** (`alpha = 1`). It does not yet support non-linear Hamiltonian cost currents or general `alpha != 1` cases.

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
    "Vertices List" -> {1, 2, 3},
    "Adjacency Matrix" -> {{0, 1, 0}, {0, 0, 1}, {0, 0, 0}},
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

Example factory workflow:

```mathematica
f = GetExampleScenario[12];
s = f[{{1, 100}}, {{4, 0}}, {}, 1, 0, Function[z, -1/z]];
```

## Active API (current phase)

Scenario kernel:
- `scenario`, `scenarioQ`
- `makeScenario`, `validateScenario`, `completeScenario`
- `ScenarioData`

Example scenario factories:
- `GetExampleScenario`
- `GridScenario`, `CycleScenario`, `GraphScenario`, `AMScenario`

Unknowns/system kernels:
- `unknowns`, `unknownsQ`, `UnknownsData`, `makeUnknowns`
- `mfgSystem`, `mfgSystemQ`, `SystemData`, `makeSystem`

Core symbolic primitives:
- `j`, `u`, `z`, `alpha`, `Cost`

## Intentionally unavailable in this phase

Not loaded by current `MFGraphs.wl` bootstrap:
- `ScenarioByKey`, `GetExampleData`
- `DataToEquations`, `CriticalCongestionSolver`, `SolveMFG`
- legacy/extended solver modules and solver benchmarking paths

## Running tests

From repository root:

```bash
wolframscript -file Scripts/RunTests.wls fast
wolframscript -file MFGraphs/Tests/scenario-kernel.mt
wolframscript -file MFGraphs/Tests/make-unknowns.mt
```

Current active runner suites (`Scripts/RunTests.wls`):
- `fast`: `scenario-kernel.mt`, `make-unknowns.mt`, `reduce-system.mt`, `scenario-consistency.mt`, `graphics.mt`
- `all`: alias for `fast`
- `archive`: archived compatibility/legacy suites (explicit use only)
- `full`: `fast + archive`

## Repository structure

```text
MFGraphs/
  MFGraphs.wl
  Scenario.wl
  Unknowns.wl
  System.wl
  Examples/
    ExampleScenarios.wl
    archive/
      ExamplesData.wl
  archive/
    DataToEquations.wl
  Tests/
    scenario-kernel.mt
    make-unknowns.mt
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
