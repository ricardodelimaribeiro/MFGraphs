# MFGraphs

A Wolfram Language package for solving **Mean Field Games on networks** with congestion and switching costs.

MFGraphs converts network topology (vertices, edges, entry/exit flows, switching costs) into systems of equations, then solves for equilibrium flow distributions and value functions using symbolic and numerical methods.

## Quick Navigation

- [Installation](#installation)
- [Quick start](#quick-start)
- [Defining a network](#defining-a-network)
- [Solvers](#solvers) — CriticalCongestionSolver, NonLinearSolver, MonotoneSolver
- [Switching costs](#switching-costs)
- [Built-in examples](#built-in-examples)
- [Configuration](#configuration)
- [Package structure](#package-structure)
- [Running tests](#running-tests)
- [Plotting results](#plotting)

For developer guidance, see [CLAUDE.md](CLAUDE.md). For troubleshooting, see [TROUBLESHOOTING.md](TROUBLESHOOTING.md).

## Repository layout

| Path | Purpose |
|---|---|
| `MFGraphs/` | Package source (`.wl` files, tests, examples) |
| `Scripts/` | Test runner, benchmarks, documentation generator |
| `Results/` | Benchmark and profiling outputs (mostly gitignored) |
| `docs/` | Internal planning and validation notes |
| `research/papers/` | Reference papers |
| `repro/` | One-off debugging and reproduction scripts |
| `DNF_PERFORMANCE_HISTORY.md` | Auto-appended by `CompareDNF.wls` — do not edit manually |
| `PARALLEL_PERFORMANCE_HISTORY.md` | Auto-appended by `BenchmarkSuite.wls` — do not edit manually |

## Installation

Clone this repository and load the package in Mathematica:

```mathematica
(* Option 1: If the package is on your $Path *)
Needs["MFGraphs`"]

(* Option 2: Direct load from a specific location *)
Get["/path/to/MFGraphs/MFGraphs/MFGraphs.wl"]
```

Requires Mathematica 12.0 or later.

## Documentation

Full API documentation is available in [API_REFERENCE.md](API_REFERENCE.md). This document is automatically generated from the package's `::usage` metadata.

## Quick start

```mathematica
(* Load the package *)
<< MFGraphs`

(* Pick a built-in test case: a 4-vertex "attraction" network *)
Data = GetExampleData[12] /. {I1 -> 100, U1 -> 0};

(* Convert network data to equations *)
d2e = DataToEquations[Data];

(* Solve the critical congestion case (all flows start at zero) *)
result = CriticalCongestionSolver[d2e];

(* The solution is stored under "AssoCritical" *)
result["AssoCritical"]
```

Stationary solver results are standardized by default. Each solver returns an association
with keys such as `"Solver"`, `"ResultKind"`, `"Feasibility"`, `"Message"`,
`"Solution"`, `"ComparableFlowVector"`, `"KirchhoffResidual"`, and the solver-specific
payload key when available: `"AssoCritical"`, `"AssoNonCritical"`, or `"AssoMonotone"`.

## Defining a network

A network is an `Association` with five required keys:

```mathematica
Data = <|
  "Vertices List"                  -> {1, 2, 3, 4},
  "Adjacency Matrix"               -> {{0,1,0,0}, {0,0,1,0}, {0,0,0,1}, {0,0,0,0}},
  "Entrance Vertices and Flows"    -> {{1, 200}},
  "Exit Vertices and Terminal Costs" -> {{4, 0}},
  "Switching Costs"                -> {{1,2,3, 5}, {3,2,1, 5}}
|>;
```

| Key | Description |
|---|---|
| `"Vertices List"` | List of vertex labels |
| `"Adjacency Matrix"` | Square matrix; entry `(i,j) = 1` means an edge from vertex `i` to vertex `j` |
| `"Entrance Vertices and Flows"` | `{{vertex, flow}, ...}` — where agents enter and how many |
| `"Exit Vertices and Terminal Costs"` | `{{vertex, cost}, ...}` — where agents exit and the terminal cost |
| `"Switching Costs"` | `{{from, at, to, cost}, ...}` — cost of switching from one edge to another at a vertex. Use `{}` for no switching costs |

Symbolic values (e.g., `I1`, `U1`, `S1`) can be used and substituted later with `/.` rules.

## Solvers

### Critical congestion solver

Solves the special case where all edge flows are zero. This is typically the fastest solver and serves as the starting point for the non-linear solver.

```mathematica
Data = GetExampleData[7] /. {I1 -> 50, U1 -> 0, U2 -> 0};
d2e = DataToEquations[Data];
result = CriticalCongestionSolver[d2e];
result["AssoCritical"]
```

### Non-linear (general congestion) solver

Iteratively solves the full non-linear problem using fixed-point iteration. Uses the critical congestion solution as an initial guess.

```mathematica
(* Define the potential function V and parameters before solving *)
V = Function[{x, edge}, 0.5 Sin[2 Pi (x + 1/4)]^2];

Data = GetExampleData[12] /. {I1 -> 100, U1 -> 0};
d2e = DataToEquations[Data];

(* Solve iteratively (default: up to 15 iterations) *)
result = NonLinearSolver[d2e];

(* Check the solution *)
IsNonLinearSolution[result]

(* Access the solution *)
result["AssoNonCritical"]
```

Options:
- `"MaxIterations"` (default `15`) — maximum number of fixed-point iterations

```mathematica
result = NonLinearSolver[d2e, "MaxIterations" -> 30];
```

### Monotone operator solver

An ODE-based solver using gradient flow on the Kirchhoff matrix.

```mathematica
Data = GetExampleData[3] /. {I1 -> 80, U1 -> 0};
solution = MonotoneSolverFromData[Data]
```

Standardized result envelope:

```mathematica
result = MonotoneSolverFromData[Data];
result["AssoMonotone"]
result["Convergence"]
```

Options:
- `"ResidualTolerance"` (default `10^-6`) — convergence threshold on the reduced monotone flow residual
- `"MaxTime"` (default `100`) — maximum pseudo-time horizon for the ODE integration
- `"MaxSteps"` (default `5000`) — maximum `NDSolveValue` step budget
- `"UseCachedProjection"` (default `True`) — reuse the projected linear solve when the state changes only slightly

```mathematica
solution = MonotoneSolverFromData[Data, "ResidualTolerance" -> 10^-6, "MaxTime" -> 20, "MaxSteps" -> 2000]
```

## Running tests

Use the suite runner for the common regression sets:

```bash
wolframscript -file Scripts/RunTests.wls fast
wolframscript -file Scripts/RunTests.wls slow
wolframscript -file Scripts/RunTests.wls all
```

Tests are currently intended to run locally via `wolframscript` commands above.

## Plotting

After solving with the non-linear solver, plot mass densities and value functions on each edge:

```mathematica
result = NonLinearSolver[d2e];

(* Plot mass density M on a specific edge *)
PlotMassDensity[result, "AssoNonCritical", {1, 2}]

(* Plot all mass densities *)
PlotMassDensities[result, "AssoNonCritical"]

(* Plot value function U on a specific edge *)
PlotValueFunction[result, "AssoNonCritical", {1, 2}]

(* Plot all value functions *)
PlotValueFunctions[result, "AssoNonCritical"]
```

## Switching costs

Switching costs model the price agents pay when transitioning between edges at a vertex. Each entry `{from, at, to, cost}` means: an agent on edge `(from, at)` switching to edge `(at, to)` pays `cost`.

```mathematica
(* Y-network with switching costs *)
Data = GetExampleData[8] /. {I1 -> 100, U1 -> 0, U2 -> 0,
                     S1 -> 2, S2 -> 3, S3 -> 2, S4 -> 1, S5 -> 3, S6 -> 1};
d2e = DataToEquations[Data];

(* Verify switching costs satisfy triangle inequality *)
IsSwitchingCostConsistent[Normal @ d2e["SwitchingCosts"]]
(* True *)
```

## Built-in examples

Use `GetExampleData[key]` to load predefined test cases:

| Key | Network |
|---|---|
| `1`–`6` | Linear chains (1 to 9 edges) |
| `7`, `8` | Y-network, 1 entrance, 2 exits (without/with switching) |
| `9`, `10` | Y-network, 2 entrances, 1 exit (without/with switching) |
| `11`, `12` | Attraction problem (with/without switching) |
| `13` | Attraction without middle edge (numerical) |
| `14` | Triangle with switching |
| `15` | Y-network, 3 vertices, 2 entrances, with switching |
| `16` | Invalid switching-cost example |
| `17`, `18` | 2-edge variants |
| `19` | Y-network with 2 entrances and switching |
| `20`–`23` | Larger multi-entrance, multi-exit networks |
| `27` | Cycle (2 vertices, undirected) |
| `104` | Triangle with 2 entrances and 3 exits |
| `105` | Numeric alias for chain with two exits |
| `"triangle with two exits"` | Triangle with 1 entrance and 2 exits |
| `"chain with two exits"` | Chain with 1 entrance and 2 exits |
| `"Braess split"` | Braess paradox — split variant |
| `"Braess congest"` | Braess paradox — congestion variant |
| `"New Braess"` | Braess with edge-dependent costs |
| `"Big Braess split"` | Extended Braess — split variant |
| `"Big Braess congest"` | Extended Braess — congestion variant |
| `"Jamaratv9"` | Jamarat pilgrimage network (9 vertices) |
| `"HRF Scenario 1"` | HRF reference scenario |
| `"Paper example"` | 4-vertex network with 2 entrances, 2 exits, and switching |
| `"Inconsistent Y shortcut"` | Y-network with deliberately inconsistent switching costs |
| `"Inconsistent attraction shortcut"` | Attraction network with a deliberately inconsistent shortcut turn |
| `"Grid0303"` | 3 x 3 grid graph |
| `"Grid0404"` | 4 x 4 grid graph |
| `"Grid0505"` | 5 x 5 grid graph |
| `"Grid0707"` | 7 x 7 grid graph |
| `"Grid0710"` | 7 x 10 grid graph |
| `"Grid1010"` | 10 x 10 grid graph |
| `"Grid1020"` | 10 x 20 grid graph |

## Configuration

### Verbose output

Progress and timing messages are printed by default. Suppress them with:

```mathematica
$MFGraphsVerbose = False;
```

### Hamiltonian parameters

The non-linear solver uses a Hamiltonian framework. Override these before calling `NonLinearSolver`:

```mathematica
(* Congestion exponent (default: 1) *)
alpha[edge_] := 1

(* Interaction potential (default: -1/m^2) *)
g[m_, edge_] := -1/m^2

(* Potential function (must be defined by the user for non-linear solving) *)
V = Function[{x, edge}, 0.5 Sin[2 Pi (x + 1/4)]^2];
```

## Package structure

```
MFGraphs/
  MFGraphs.wl             Package loader, SolveMFG unified entrypoint
  DNFReduce.wl            Boolean algebra (disjunctive normal form reduction)
  DataToEquations.wl      Network topology → equation converter
  Solvers.wl              Critical-congestion solver suite (extracted Phase 3)
  NonLinearSolver.wl      Iterative non-linear solver and Hamiltonian framework
  Monotone.wl             Monotone operator (ODE-based) solver
  TimeDependentSolver.wl  Time-dependent MFG solver
  Graphics.wl             Public visualization helpers (NetworkGraphPlot, SolutionFlowPlot, ExitFlowPlot)
  Examples/
    ExamplesData.wl        34 built-in test cases via GetExampleData[key]
  Tests/
    *.mt                   MUnit test files
  Kernel/
    init.m                 Paclet initialization
Scripts/
  RunTests.wls             Test suite runner
  BenchmarkSuite.wls       Performance benchmarking tool
  CompareDNF.wls           DNF solver before/after comparison
  GenerateDocs.wls         Automated API documentation generator
```

### Maintenance

To regenerate the API documentation:
```bash
wolframscript -file Scripts/GenerateDocs.wls
```
```

## Examples

All examples are reproducible by running the code in this README using the `GetExampleData[key]` function, which provides 34 built-in test cases. See the [Built-in examples](#built-in-examples) table for available keys.

## End-to-end example

A complete workflow for a Braess-paradox network with congestion:

```mathematica
<< MFGraphs`

(* Load the Braess congestion example *)
Data = GetExampleData["Braess congest"];

(* Convert to equations *)
d2e = DataToEquations[Data];

(* View the constructed graph *)
d2e["auxiliaryGraph"]

(* Solve the critical congestion case *)
critical = CriticalCongestionSolver[d2e];
critical["AssoCritical"]

(* Define potential and solve the non-linear case *)
V = Function[{x, edge}, 0];

nonlinear = NonLinearSolver[critical];
IsNonLinearSolution[nonlinear]

(* Plot the results *)
GraphicsGrid[Partition[PlotMassDensities[nonlinear, "AssoNonCritical"], 3]]
GraphicsGrid[Partition[PlotValueFunctions[nonlinear, "AssoNonCritical"], 3]]
```

## License

This project is part of ongoing research. Please contact the authors before using in publications.
