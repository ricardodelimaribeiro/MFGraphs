# MFGraphs

A Wolfram Language package for solving **Mean Field Games on networks** with congestion and switching costs.

MFGraphs converts network topology (vertices, edges, entry/exit flows, switching costs) into systems of equations, then solves for equilibrium flow distributions and value functions using symbolic and numerical methods.

The model keeps its Hamiltonian structure. The active runtime solver surface in this repository targets the critical specialization with edgewise \(\alpha=1\).

## Quick Navigation

- [Installation](#installation)
- [Quick start](#quick-start)
- [Defining a network](#defining-a-network)
- [Solvers](#solvers)
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
| `docs/history/DNF_PERFORMANCE_HISTORY.md` | Auto-appended by `CompareDNF.wls` — do not edit manually |
| `docs/history/PARALLEL_PERFORMANCE_HISTORY.md` | Auto-appended by `BenchmarkSuite.wls` — do not edit manually |

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

Stationary solver results are standardized. Each solver returns an association
with keys such as `"Solver"`, `"ResultKind"`, `"Feasibility"`, `"Message"`,
`"Solution"`, `"ComparableFlowVector"`, and `"KirchhoffResidual"`.

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

Solves the critical congestion specialization of the Hamiltonian model ($\alpha=1$ on every edge). This solver uses a combination of symbolic DNF reduction and numeric iterative backends (Fictitious Play) to find global equilibria.

```mathematica
Data = GetExampleData[7] /. {I1 -> 50, U1 -> 0, U2 -> 0};
d2e = DataToEquations[Data];
result = CriticalCongestionSolver[d2e];
result["AssoCritical"]
```

## Running tests

Use the suite runner for the common regression sets:

```bash
wolframscript -file Scripts/RunTests.wls fast
wolframscript -file Scripts/RunTests.wls slow
```

## Plotting

Visualize the network and the resulting flows:

```mathematica
result = CriticalCongestionSolver[d2e];

(* Plot net flows on the network graph *)
SolutionFlowPlot[d2e, result["Solution"]]

(* Plot exit flow distribution *)
ExitFlowPlot[result["ExitFlows"]]
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

Use `GetExampleData[key]` to load predefined test cases (34 available). Common keys include:

| Key | Network |
|---|---|
| `7`, `8` | Y-network, 1 entrance, 2 exits (without/with switching) |
| `11`, `12` | Attraction problem (with/without switching) |
| `"Braess congest"` | Braess paradox — congestion variant |
| `"Jamaratv9"` | Jamarat pilgrimage network (9 vertices) |
| `"Paper example"` | 4-vertex network with 2 entrances, 2 exits, and switching |

## Configuration

### Verbose output

Progress and timing messages are printed by default. Suppress them with:

```mathematica
$MFGraphsVerbose = False;
```

## Package structure

```
MFGraphs/
  MFGraphs.wl             Package loader, SolveMFG unified entrypoint
  DNFReduce.wl            Boolean algebra (disjunctive normal form reduction)
  DataToEquations.wl      Network topology → equation converter
  Solvers.wl              Critical-congestion solver suite
  Graphics.wl             Public visualization helpers (NetworkGraphPlot, SolutionFlowPlot, ExitFlowPlot)
  Examples/
    ExamplesData.wl        Built-in test cases via GetExampleData[key]
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

## Archived non-critical solvers

Non-critical solver families are archived from active runtime paths, but recoverable from tags:

- `archive/non-critical-solvers`
- `archive/non-critical-solvers-full`

Useful retrieval commands:

```bash
git show archive/non-critical-solvers:MFGraphs/NonLinearSolver.wl
git checkout archive/non-critical-solvers -- MFGraphs/NonLinearSolver.wl
git show archive/non-critical-solvers-full:MFGraphs/Monotone.wl
git checkout archive/non-critical-solvers-full -- MFGraphs/Monotone.wl
```

## License

This project is part of ongoing research. Please contact the authors before using in publications.
