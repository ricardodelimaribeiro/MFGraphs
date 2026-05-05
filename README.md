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

### Signed edge traversal cost

Physical network edges are undirected, but the symbolic system keeps one flow variable per direction. For an edge `{a, b}`, the oriented net flow is:

```mathematica
J = j[a, b] - j[b, a]
```

`buildHamiltonianData` stores the Hamiltonian right-hand side metadata as:

```mathematica
J - Sign[J] J^alpha
```

For critical congestion (`Alpha == 1`), the `Sign[J]` factor makes this orientation-aware:

```mathematica
J - Sign[J] J
```

and the congestion-cost part is equivalent to:

```mathematica
J Sign[J] == Abs[J]
```

Thus expressions such as:

```mathematica
-(j[1, 5] - j[5, 1]) Sign[j[1, 5] - j[5, 1]]
```

are the signed form of the traversal cost on the undirected edge `{1, 5}`. The sign convention lets one formula handle either travel direction along the edge. Current critical solvers solve the critical `EqGeneral` constraints directly; `Nrhs` and `CostCurrents` remain stored Hamiltonian metadata for diagnostics and future non-critical/numeric paths.

### Solver diagnostics proof of concept

Legacy-inspired diagnostics and alternate critical solvers are opt-in. They do not change `solveScenario`, `dnfReduceSystem`, or raw solver return shapes.

```mathematica
Needs["MFGraphs`"];

s = gridScenario[{3}, {{1, 120}}, {{2, 0}, {3, 10}}];
sys = makeSystem[s];

defaultSol = solveScenario[s];
report = solutionReport[sys, defaultSol];

report["ResultKind"]
report["ValidationReport"]["Valid"]
report["KirchhoffResidual"]
report["TransitionFlowDiagnostics"]

flowFirstSol = flowFirstCriticalSystem[sys];
isValidSystemSolution[sys, flowFirstSol]
```

For residual outputs, branch costs can be ranked without rewriting the solution:

```mathematica
branchScenario = getExampleScenario[12, {{1, 100}}, {{4, 0}}];
branchSys = makeSystem[branchScenario];
branchedSol = optimizedDNFReduceSystem[branchSys];

branchReport = solutionBranchCostReport[branchSys, branchedSol];

branchReport["BranchCount"]
KeyTake[#, {"BranchIndex", "TotalObjective"}] & /@ branchReport["BestBranches"]
```

`directCriticalSystem[sys]` is stricter: it only runs on critical systems where every switching cost is numeric zero and all exit costs are equal numeric values. Otherwise it returns `Failure["directCriticalSystem", <|"Reason" -> ...|>]`.

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
- `solutionReport`, `solutionBranchCostReport`, `directCriticalSystem`, `flowFirstCriticalSystem`

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
