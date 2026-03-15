# MFGraphs

A Wolfram Language package for solving **Mean Field Games on networks** with congestion and switching costs.

MFGraphs converts network topology (vertices, edges, entry/exit flows, switching costs) into systems of equations, then solves for equilibrium flow distributions and value functions using symbolic and numerical methods.

## Installation

Clone this repository and load the package in Mathematica:

```mathematica
(* Option 1: If the package is on your $Path *)
Needs["MFGraphs`"]

(* Option 2: Direct load from a specific location *)
Get["/path/to/MFGraphs/MFGraphs/MFGraphs.m"]
```

Requires Mathematica 12.0 or later.

## Quick start

```mathematica
(* Load the package *)
<< MFGraphs`

(* Pick a built-in test case: a 4-vertex "attraction" network *)
Data = DataG[12] /. {I1 -> 100, U1 -> 0};

(* Convert network data to equations *)
d2e = Data2Equations[Data];

(* Solve the critical congestion case (all flows start at zero) *)
result = CriticalCongestionSolver[d2e];

(* The solution is stored under "AssoCritical" *)
result["AssoCritical"]
```

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
Data = DataG[7] /. {I1 -> 50, U1 -> 0, U2 -> 0};
d2e = Data2Equations[Data];
result = CriticalCongestionSolver[d2e];
result["AssoCritical"]
```

### Non-linear (general congestion) solver

Iteratively solves the full non-linear problem using fixed-point iteration. Uses the critical congestion solution as an initial guess.

```mathematica
(* Define the potential function V and parameters before solving *)
V = Function[{x, edge}, 0.5 Sin[2 Pi (x + 1/4)]^2];

Data = DataG[12] /. {I1 -> 100, U1 -> 0};
d2e = Data2Equations[Data];

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
Data = DataG[3] /. {I1 -> 80, U1 -> 0};
solution = MonotoneSolverFromData[Data]
```

Options:
- `"TimeSteps"` (default `100`) — number of ODE integration steps

```mathematica
solution = MonotoneSolverFromData[Data, "TimeSteps" -> 200]
```

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
Data = DataG[8] /. {I1 -> 100, U1 -> 0, U2 -> 0,
                     S1 -> 2, S2 -> 3, S3 -> 2, S4 -> 1, S5 -> 3, S6 -> 1};
d2e = Data2Equations[Data];

(* Verify switching costs satisfy triangle inequality *)
IsSwitchingCostConsistent[Normal @ d2e["SwitchingCosts"]]
(* True *)
```

## Built-in examples

Use `DataG[key]` to load predefined test cases:

| Key | Network |
|---|---|
| `1`–`6` | Linear chains (1 to 9 edges) |
| `7`, `8` | Y-network, 1 entrance, 2 exits (without/with switching) |
| `9`, `10` | Y-network, 2 entrances, 1 exit (without/with switching) |
| `11`, `12` | Attraction problem (with/without switching) |
| `13` | Attraction without middle edge (numerical) |
| `14` | Triangle with switching |
| `15` | Y-network, 3 vertices, 2 entrances, with switching |
| `17`, `18` | 2-edge variants |
| `19` | Y-network with 2 entrances and switching |
| `20`–`23` | Larger multi-entrance, multi-exit networks |
| `27` | Cycle (2 vertices, undirected) |
| `104` | Triangle with 2 entrances and 3 exits |
| `"Braess split"` | Braess paradox — split variant |
| `"Braess congest"` | Braess paradox — congestion variant |
| `"New Braess"` | Braess with edge-dependent costs |
| `"Big Braess split"` | Extended Braess — split variant |
| `"Big Braess congest"` | Extended Braess — congestion variant |
| `"Jamaratv9"` | Jamarat pilgrimage network (9 vertices) |
| `"Paper example"` | 4-vertex network with 2 entrances, 2 exits, and switching |
| `"Grid0303"` | 3 x 3 grid graph |
| `"Grid1020"` | 10 x 20 grid graph |
| `"triangle with two exits"` | Triangle with 1 entrance and 2 exits |

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
  MFGraphs.m              Package loader, verbose flag, MFGPrint helpers
  DNFReduce.m             Boolean algebra (disjunctive normal form reduction)
  D2E2.m                  Data-to-equations converter and critical congestion solver
  NonLinearSolver.m       Iterative non-linear solver and Hamiltonian framework
  Monotone.m              Monotone operator (ODE-based) solver
  Examples/
    ExamplesData.m         Built-in test cases and DataG function
  Kernel/
    init.m                 Paclet initialization
```

## Examples

All examples are reproducible by running the code in this README using the `DataG[key]` function, which provides 34 built-in test cases. See the [Built-in examples](#built-in-examples) table for available keys.

## End-to-end example

A complete workflow for a Braess-paradox network with congestion:

```mathematica
<< MFGraphs`

(* Load the Braess congestion example *)
Data = DataG["Braess congest"];

(* Convert to equations *)
d2e = Data2Equations[Data];

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
