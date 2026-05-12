# MFGraphs

**MFGraphs** is a Wolfram Language package for Mean Field Games on networks.

**New to the package?** See [TOUR.md](TOUR.md) for a guided reading order.

The repository is currently in a **core scenario-kernel phase**. The active package surface is centered on:
- typed scenario construction
- example scenario factories
- exact symbolic unknown bundle construction
- structural system construction
- critical-congestion symbolic solving
- visualization helpers

## Mean Field Games on Networks: Model Overview

MFGraphs solves for a stationary Mean Field Game (MFG) equilibrium on a network graph. The model describes the strategic behavior of a large population of rational agents moving through a graph.

**Key Components:**
1. **Mass Conservation:** A stationary distribution of agents satisfying Kirchhoff's laws (flow in = flow out) at all **original vertices** of the graph. Boundary conditions (entries/exits) are handled by auxiliary vertices.
2. **Congestion Model:** The model employs an explicit Hamiltonian:
   $$H(x,p,z) = \frac{|p|^2}{2 z^\alpha} + V(x) - g(z)$$
   where $p$ is the momentum, $z$ is the density, $V(x)$ is the potential, and $g(z)$ represents the congestion cost.
3. **Fundamental Equations:**
   - **Hamilton-Jacobi-Bellman (HJB) Equation:** $H(x, u_x, m) = 0$, where $u$ is the value function and $m$ is mass.
   - **Transport Equation:** $(-m D_p H(x, u_x, m))_x = 0$.
4. **Edge Flow Relation:** Based on the model Hamiltonian, the flow is given by $j = -m \frac{u_x}{m^\alpha}$. Integrating over a canonical edge segment $[0, 1]$ for the **critical congestion** case ($\alpha = 1$) yields the constant flow relation:
   $$j = u(1) - u(0)$$
5. **Equilibrium (Complementarity):** Flow only exists on edges that are part of an optimal path, resulting in a Linear Complementarity Problem (LCP) for the critical regime.

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

The solver layer combines package-defined MFGraphs reductions with Wolfram
backends:

| Solver | Method |
|---|---|
| `reduceSystem` | Wolfram `Reduce` over `Reals` (CAD / real quantifier elimination for polynomial real systems). |
| `dnfReduceSystem` | Package recursive DNF reducer: equality substitution and disjunction expansion. |
| `optimizedDNFReduceSystem` | Package branch-state DNF reduction for small residual systems, falling back to `dnfReduce`. |
| `activeSetReduceSystem` | Package active-set enumeration for complementarity alternatives, falling back to `dnfReduce`. |
| `booleanReduceSystem` | Wolfram `BooleanConvert[..., "DNF"]`, then `Reduce` per disjunct. |
| `booleanMinimizeSystem` | Wolfram `BooleanMinimize[..., "DNF"]`, then `Reduce` per disjunct. This is exact Boolean minimal-DNF/SOP minimization; it is analogous to classical Quine-McCluskey/Petrick-style two-level minimization, but Wolfram does not document `BooleanMinimize` as specifically using QMC or Petrick internally. |
| `booleanMinimizeReduceSystem` | Package arm pruning and component decomposition, then `BooleanMinimize` and `Reduce` per component/disjunct. |
| `findInstanceSystem` | Wolfram `FindInstance` over `Reals`. |
| `directCriticalSystem` | Package graph-distance pruning heuristic, then `FindInstance`. |
| `flowFirstCriticalSystem` | Package linear least-squares candidate via `LeastSquares`, then `FindInstance` fallback. |

See [Solver Methods and Package Novelty](docs/solvers/methods-and-novelty.md)
for the full method map and the distinction between Wolfram backend algorithms
and MFGraphs-specific contributions.

## Package novelty

The package novelty is the MFG-specific workflow around those solver backends:
typed scenario construction, auxiliary state-space graph augmentation, exact
logical-algebraic complementarity encoding, equality preprocessing,
custom DNF/active-set reductions, hybrid Boolean/real branch solving, solver
diagnostics, and paired raw/rich network visualization. The package should not
claim novelty for Wolfram built-ins such as `Reduce`, `FindInstance`,
`BooleanMinimize`, `RowReduce`, or `LeastSquares`.

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
- Solver methods and novelty: [docs/solvers/methods-and-novelty.md](docs/solvers/methods-and-novelty.md)
- Troubleshooting: [TROUBLESHOOTING.md](TROUBLESHOOTING.md)

## License

This project is part of ongoing research. Please contact the authors before using it in publications.
