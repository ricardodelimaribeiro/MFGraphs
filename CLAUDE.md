# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**MFGraphs** is a Wolfram Language paclet (`PacletInfo.m`, v0.0.2) for solving stationary **Mean Field Games on networks**. A large population of rational agents move through a graph under Hamiltonian dynamics; mass conservation (Kirchhoff laws) and complementarity constraints yield a Linear Complementarity Problem that the package builds and reduces symbolically.

The repository is currently oriented toward an in-progress paper on the **critical congestion solver** (the Alpha = 1 regime). The active package surface — scenario construction, symbolic unknowns, structural system assembly, and critical-congestion symbolic solvers — is exactly what supports that paper. Anything outside that scope lives in `MFGraphs/archive/` and is **not** loaded by `Needs["MFGraphs`"]`.

## Environment

- Mathematica / Wolfram Language **12.0+**
- `wolframscript` on PATH; local app at `/Applications/Wolfram.app`
- Quick check: `wolframscript -version`

## Commands

```bash
# Load the package interactively
wolframscript -e 'Needs["MFGraphs`"]'

# Run the active test suite
wolframscript -file Scripts/RunTests.wls fast

# Run a single test file
wolframscript -file Scripts/RunSingleTest.wls MFGraphs/Tests/scenario-kernel.mt

# Benchmark the system solver (required before/after for solver-sensitive PRs)
wolframscript -file Scripts/BenchmarkSystemSolver.wls --tag "short-description" --timeout 60

# Profile the scenario kernel (non-mutating stage timings)
wolframscript -file Scripts/ProfileScenarioKernel.wls --case example-12 --timeout 10

# Profile the DNF reducer
wolframscript -file Scripts/ProfileDNFReduce.wls --case all --order all --timeout 60

# Regenerate API_REFERENCE.md from ::usage strings (do not edit it by hand)
wolframscript -file Scripts/GenerateDocs.wls
```

## Architecture

**Loading order and dependencies.** `Needs["MFGraphs`"]` loads each sub-package as its own flat context, then `EndPackage` flattens everything onto `$ContextPath` so public symbols are unqualified. `MFGraphs/archive/` is intentionally excluded.

```
primitives  →  utilities  →  scenarioTools  →  examples
                                            →  unknownsTools  →  systemTools  →  Tawaf
                                                                              →  solversTools  →  orchestrationTools
                                                                                              →  graphicsTools
```

**Composition pipeline.** Every workflow is the same three-step chain, optionally followed by orchestration:

```
makeScenario  →  makeSymbolicUnknowns  →  makeSystem  →  solveScenario
```

- `makeScenario` normalizes input (Graph → adjacency), validates, fills defaults, and caches `buildAuxiliaryTopology` output inside the returned `scenario[<|...|>]`.
- `makeSymbolicUnknowns` reads the cached topology and emits the symbolic variable families `j`, `u`, `z`.
- `makeSystem` runs `buildBoundaryData`, `buildFlowData`, `buildComplementarityData`, `buildHamiltonianData` and merges their typed records into one `mfgSystem`.
- `solveScenario` (in `orchestrationTools`) wraps the chain end-to-end.

**Typed wrapper pattern.** Every domain object is `Head[Association]` plus a predicate and an accessor. `typeData[obj, key]` returns `Missing["KeyAbsent", key]` if the key is absent.

| Head | Predicate | Accessor |
|---|---|---|
| `scenario` | `scenarioQ` | `scenarioData[s]` / `scenarioData[s, key]` |
| `symbolicUnknowns` | `symbolicUnknownsQ` | `symbolicUnknownsData[u]` / `symbolicUnknownsData[u, key]` |
| `mfgSystem` | `mfgSystemQ` | `systemData[sys]` / `systemData[sys, key]` |
| `mfgBoundaryData` / `mfgFlowData` / `mfgComplementarityData` / `mfgHamiltonianData` | typed sub-records nested inside `mfgSystem` | use `systemData[sys, key]` or `systemDataFlatten[sys]` |

**Naming conventions.** Enforced throughout the package:

| Prefix | Role |
|---|---|
| `make*` | Eager constructor: validate + complete + wrap in typed head |
| `build*` | Structural builder: compute topology / matrices / equations |
| `*Q` | Structural predicate |
| `validate*` | Spec validator (returns True/False with diagnostics) |

## Performance Patterns

- **Topology caching.** `buildAuxiliaryTopology` runs once inside `makeScenario`; downstream builders read from the cached `scenario[]` rather than recomputing.
- **DNF-first reduction.** `dnfReduceSystem` / `optimizedDNFReduceSystem` factor the complementarity system into a disjunction of branches. The optimized variant prunes infeasible branches early; the plain variant is the reference.
- **Active-set reducer.** `activeSetReduceSystem` is the alternative path used when DNF-first is too expensive on a given topology.
- **Verbose gating.** Use `mfgPrint[...]` (not `Print`) for debug output — it is silent unless `$MFGraphsVerbose = True`.

## Test Suites

Active `.mt` files under `MFGraphs/Tests/`:

| Suite | Covers |
|---|---|
| `scenario-kernel.mt` | `makeScenario` validation, defaults, topology caching |
| `symbolic-unknowns.mt` | `makeSymbolicUnknowns` variable bundle generation |
| `reduce-system.mt` | system reduction end-to-end |
| `scenario-consistency.mt` | cross-stage invariants (scenario ↔ unknowns ↔ system) |
| `dnf-reducer.mt` | DNF reducer correctness and edge cases |
| `boolean-minimize.mt` | `booleanMinimizeSystem` / `booleanMinimizeReduceSystem` |
| `orchestration.mt` | `solveScenario` and full pipeline |
| `graphicsTools.mt` | `rawNetworkPlot`, `richNetworkPlot` |
| `tawaf.mt` | unrolled-circumambulation scenario builder |

## Conventions

- Avoid single-letter System symbol shadows (`K`, `D`, `C`, ...) — use descriptive names so the `MFGraphs`` context does not collide with built-ins.
- Prefer typed constructors and accessors over building raw associations by hand.
- **Solver-sensitive PRs must include before/after benchmark output** produced with `BenchmarkSystemSolver.wls --tag`. Append results to `BENCHMARKS.md`.
- Do not commit to `master`; use feature branches and PRs.
- `API_REFERENCE.md` is generated from `::usage` strings — edit the source, then run `Scripts/GenerateDocs.wls`.

## Documentation

- `README.md` — model overview and mathematical setup
- `CONTRIBUTING.md` — PR conventions, commit style, branch policy
- `API_REFERENCE.md` — auto-generated public API reference (do not edit)
- `BENCHMARKS.md` — solver benchmark policy and historical tagged results
- `TROUBLESHOOTING.md` — package load failures, scenario validation errors, solver timeouts

## Project Skills

This repo ships three Claude Code skills: `/test` (run the fast suite and report pass/fail), `/scenario` (walk `makeScenario → makeSystem → ReduceSystem` with consistency checks), and `/ship` (auto-branch, test, commit, push, open PR). Prefer them over reinventing the equivalent shell pipelines.
