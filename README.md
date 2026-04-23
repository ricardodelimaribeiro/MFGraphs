# MFGraphs

**MFGraphs** is a Wolfram Language package for studying **Mean Field Games on networks** with congestion and switching costs. The current repository is organized around a **critical-congestion runtime**, a growing **scenario kernel**, a set of plotting and helper utilities, and an active documentation/planning effort aimed at a **scenario-first public API**.

The most important practical point is that the repository is **not just the `v0.4.0-critical-only` milestone anymore**. That tag remains a useful historical checkpoint, but current `master` contains additional modules, tests, planning documents, and helper layers beyond that milestone. The README therefore describes the **current repository state**, while also making clear which parts of the surface are already stable and which are still being normalized.

## Quick navigation

- [Repository status](#repository-status)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Canonical workflows](#canonical-workflows)
- [Core package surface](#core-package-surface)
- [Scenario and unknowns layers](#scenario-and-unknowns-layers)
- [Built-in examples](#built-in-examples)
- [Running tests](#running-tests)
- [Repository structure](#repository-structure)
- [Documentation and planning artifacts](#documentation-and-planning-artifacts)
- [Historical tags](#historical-tags)
- [License](#license)

For contributor workflow guidance, see [CLAUDE.md](CLAUDE.md). For troubleshooting notes, see [TROUBLESHOOTING.md](TROUBLESHOOTING.md).

## Repository status

The repository currently combines a **stable critical solver path** with an **ongoing scenario-first API normalization effort**. In practical terms, the safest current path is still **raw example or model data -> `DataToEquations` -> `CriticalCongestionSolver` or `SolveMFG`**. The scenario and unknowns layers are real, loaded modules, but they should still be read as **active transition infrastructure** rather than as a fully frozen long-term API.

| Area | Current status |
|---|---|
| Main runtime focus | **Critical-congestion** solver path |
| High-level entrypoint | `SolveMFG` |
| Compilation entrypoint | `DataToEquations` |
| Main solver | `CriticalCongestionSolver` |
| Typed workflow layers | `Scenario.wl`, `Unknowns.wl` are loaded and usable, but still under active API normalization |
| Plotting support | `Graphics.wl` public plotting helpers |
| Repository-only solver infrastructure | Additional files such as `FictitiousPlayBackend.wl` and `SolveMFGDispatch.wl` exist in the repo, but are not part of the default module load path described by `MFGraphs.wl` |
| Public API state | **In transition** toward a scenario-first design |
| Historical checkpoint | `v0.4.0-critical-only` remains a milestone tag, not a full description of current `master` |

## Installation

Clone the repository and load the package from Mathematica or WolframScript.

```mathematica
(* If the package directory is already on $Path *)
Needs["MFGraphs`"]

(* Or load directly from a checkout *)
Get["/path/to/MFGraphs/MFGraphs/MFGraphs.wl"]
```

The package targets **Mathematica 12.0 or later**.

## Quick start

The shortest current workflow is still **example data -> compiled equations -> critical solve**.

```mathematica
<< MFGraphs`

data = GetExampleData[12] /. {I1 -> 100, U1 -> 0};
d2e = DataToEquations[data];
result = CriticalCongestionSolver[d2e];

result["Feasibility"]
result["Solution"]
```

If you want the package-managed dispatch path, use `SolveMFG`.

```mathematica
result1 = SolveMFG[data];
result2 = SolveMFG[d2e];
```

Solver outputs are standardized associations. The exact payload varies by solver/backend, but several keys are common enough to treat as the main inspection surface.

| Common key | Meaning |
|---|---|
| `"Solver"` | Solver/backend name used for the run |
| `"ResultKind"` | High-level result classification |
| `"Feasibility"` | Feasibility verdict such as `"Feasible"` or `"Infeasible"` |
| `"Message"` | Short diagnostic summary |
| `"Solution"` | Primary solution association when available |
| `"ComparableFlowVector"` | Numeric flow representation for comparisons |
| `"KirchhoffResidual"` | Residual measure for flow conservation checks |
| `"AssoCritical"` | Critical-solver payload association for symbolic/critical workflows |

## Canonical workflows

The repository currently supports **two practical ways of working**. The first is the established raw-data path, which is still the most direct route for existing examples, regression tests, and the clearest reproducible solver behavior. The second is the newer scenario-oriented path, which is intended to become the canonical public workflow over time but should still be treated as a **transition target** rather than as the sole recommended contract.

| Workflow | Current role |
|---|---|
| Raw data -> `DataToEquations` -> solver | **Most stable current workflow** and the path best aligned with existing examples and tests |
| Scenario -> validation/completion -> compile/solve | **Emerging workflow** intended to become the long-term public shape |

At the moment, the README still shows the raw-data path first because it is the most stable and already widely exercised. However, the repository planning documents in `docs/planning/` now treat **scenario-first normalization** as the architectural direction for future public API work.

## Core package surface

The main package loader in `MFGraphs/MFGraphs.wl` currently loads a **specific default runtime set**. That loaded surface is narrower than the full set of Wolfram source files present in the repository, so it is worth distinguishing those two concepts explicitly.

| Loaded by default via `MFGraphs.wl` | Role |
|---|---|
| `Examples/ExamplesData.wl` | Built-in example and benchmark data |
| `DNFReduce.wl` | Symbolic logical reduction and solver helpers |
| `Scenario.wl` | Typed scenario kernel |
| `Unknowns.wl` | Construction helpers for symbolic unknown families |
| `System.wl` | System-kernel layer now loaded by default and used by `DataToEquations` |
| `DataToEquations.wl` | Raw network/specification data -> compiled equation association |
| `Solvers.wl` | Critical solver path and related solver utilities |
| `Graphics.wl` | Public plotting helpers |

| Present in the repository but not loaded by the default `MFGraphs.wl` sequence | Current status |
|---|---|
| `SolveMFGDispatch.wl` | Repository implementation artifact relevant to dispatch/history work, but not part of the current default load sequence |
| `FictitiousPlayBackend.wl` | Repository backend work present for numeric/backend development, but not part of the current default load sequence |

Several additional exported helpers exist today, including flow-feasibility checks, comparison builders, reduced-coordinate helpers, and backend utilities. They should not all be read as equally stable long-term API. The repository is actively distinguishing between **core workflow functions**, **advanced helpers**, **compatibility exports**, and **symbols likely to be retired or demoted** as the scenario-first design matures.

## Scenario and unknowns layers

The repository now includes a lightweight typed scenario layer together with an unknown-family construction layer. These modules are important because they represent the direction in which the package is moving: away from purely ambient symbolic manipulation and toward more structured objects with clearer validation boundaries. At the same time, the README should be explicit that this layer is **not yet the sole authoritative workflow** for the package.

```mathematica
sc = makeScenario[<|
  "Identity" -> <|"name" -> "demo"|>,
  "Model" -> data,
  "Data" -> {I1 -> 100, U1 -> 0}
|>];

scenarioQ[sc]
ScenarioData[sc, "Identity"]

unks = makeUnknowns[sc];
UnknownsData[unks, "js"]
```

| Function family | Purpose |
|---|---|
| `makeScenario`, `validateScenario`, `completeScenario`, `scenarioQ`, `ScenarioData` | Build, validate, complete, and inspect scenario objects |
| `makeUnknowns`, `unknownsQ`, `UnknownsData` | Build and inspect grouped symbolic unknown families derived from scenarios or model data |

This layer is still evolving. The active planning work assumes that future documentation should eventually present **scenario construction, validation, solve, inspect, and plot** as the main user journey, but current users should still expect the raw-data path to be the most battle-tested route.

## Built-in examples

Use `GetExampleData[key]` to retrieve predefined benchmark and teaching examples stored in `MFGraphs/Examples/ExamplesData.wl`. These remain the fastest way to exercise the package and the test suite.

| Example key | Description |
|---|---|
| `7`, `8` | Y-network stress cases without and with switching |
| `11`, `12` | Attraction examples without and with switching |
| `"Braess congest"` | Braess-type congestion example |
| `"Jamaratv9"` | Large Jamarat network variant |
| `"Paper example"` | Four-vertex example used in paper-oriented workflows |

Many examples still use symbolic parameters such as `I1`, `U1`, and `S1`. Those symbols are still present in the repository today, but part of the current API-normalization effort is to decide which of them should remain visible, which should become compatibility-only, and which should disappear from the main public workflow once scenarios carry parameter data more explicitly. Upstream test and script work has also started reducing active reliance on `GetExampleData` for routine automation, so examples should now be read primarily as **user-facing convenience data and migration material**, not as the whole future architecture.

## Running tests

For routine regression checks, use the centralized test runner from the repository root.

```bash
wolframscript -file Scripts/RunTests.wls fast
wolframscript -file Scripts/RunTests.wls slow
wolframscript -file Scripts/RunTests.wls all
wolframscript -file Scripts/RunTests.wls legacy
wolframscript -file Scripts/RunTests.wls full
```

The current repository contains an active MUnit inventory together with a newer distinction between **active scenario-first gates** and **legacy migration suites**.

| Test area | Current state |
|---|---|
| Active `.mt` files in `MFGraphs/Tests/` | **17** active test files present in the tree |
| Archived `.mt` files in `MFGraphs/Tests/archive/` | **5** archived test files plus archive README |
| Active runner gates | `fast` and `slow` now emphasize scenario-first/public-path checks |
| Legacy migration gates | `legacy-fast`, `legacy-slow`, `legacy`, and `full` keep older `GetExampleData`-tied workflows visible during migration |
| Utility runners | `Scripts/RunTests.wls`, `Scripts/RunSingleTest.wls`, `Scripts/SmokeTestExactMode.wls` |

Representative direct test runs look like this:

```bash
wolframscript -file MFGraphs/Tests/solver-contracts.mt
wolframscript -file MFGraphs/Tests/scenario-kernel.mt
wolframscript -file MFGraphs/Tests/graphics-public-api.mt
```

## Repository structure

The repository is broader than the runtime package alone. In addition to the package source, it contains benchmarking, profiling, reproduction scripts, planning notes, and generated API audit artifacts.

| Path | Purpose |
|---|---|
| `MFGraphs/` | Package source, examples, tests, kernel init, and supporting modules |
| `Scripts/` | Test runners, benchmark drivers, profiling tools, doc generation, and API-audit utilities |
| `Results/` | Benchmark and profiling outputs, mostly generated artifacts |
| `docs/` | Planning notes, validation logs, and historical documentation |
| `repro/` | Reproduction and focused verification scripts |
| `research/` | Reference papers and research support material |
| `History/` | Repository history and supporting notes |

The current runtime-oriented source tree looks like this.

```text
MFGraphs/
  MFGraphs.wl                 Main package loader and top-level usages
  DNFReduce.wl                Symbolic reduction engine
  DataToEquations.wl          Graph/specification -> compiled system
  Solvers.wl                  Critical solver path and related helpers
  SolveMFGDispatch.wl         Dispatch-related logic present in the repository
  Scenario.wl                 Typed scenario kernel
  Unknowns.wl                 Typed unknown-family construction helpers
  System.wl                   System-kernel layer extracted from DataToEquations
  Graphics.wl                 Public plotting helpers
  FictitiousPlayBackend.wl    Numeric/backend support work present in the repo
  Documentation/             Additional package-side documentation assets
  Getting started.wl         Workbook-style starter material present in the package tree
  Examples/
    ExamplesData.wl           Built-in examples
  Tests/
    *.mt                      Active MUnit regressions
    archive/                  Archived legacy tests
  Kernel/
    init.m                    Paclet initialization
```

## Documentation and planning artifacts

The repository includes both generated reference material and active planning documents. The generated API reference describes the **currently exported symbol surface**, while the planning documents describe the intended **scenario-first stabilization path**. These two views should not be conflated: exported today does **not** automatically mean stable forever.

| File or directory | Purpose |
|---|---|
| `API_REFERENCE.md` | Generated reference snapshot from `::usage` strings; useful as an export inventory, but not a curated stability guarantee |
| `Scripts/GenerateDocs.wls` | Regenerates API reference material |
| `Scripts/AuditPublicAPI.py` | Generates public-surface inventory and gap reports |
| `docs/planning/README.md` | Index for planning documents in the repository |
| `docs/planning/public_api_symbol_inventory.md` | Inventory of current exported symbols |
| `docs/planning/public_api_gap_report.md` | Documentation/test gap analysis for exported symbols |
| `docs/planning/public_api_full_surface_plan.md` | Scenario-first roadmap for API stabilization |
| `docs/planning/phase1_scenario_transition_matrix.md` | Scenario-first transition work product for API normalization |

This distinction matters. The **generated API reference** reflects what is exported today, while the **planning documents** explain why the repository is not treating every current export as equally permanent.

## Historical tags

Historical tags and archive refs remain useful for recovering earlier solver states.

| Ref | Purpose |
|---|---|
| `v0.4.0-critical-only` | Milestone tag for the earlier critical-only package-surface cleanup |
| `archive/non-critical-solvers` | Archived non-critical solver state |
| `archive/non-critical-solvers-full` | Fuller archived non-critical solver state |
| `archive/nonlinear-residual-bug` | Historical diagnostic tag related to nonlinear residual behavior |

Representative retrieval commands are shown below.

```bash
git checkout v0.4.0-critical-only
git show archive/non-critical-solvers:MFGraphs/NonLinearSolver.wl
git checkout archive/non-critical-solvers -- MFGraphs/NonLinearSolver.wl
```

## License

This project is part of ongoing research. Please contact the authors before using it in publications.
