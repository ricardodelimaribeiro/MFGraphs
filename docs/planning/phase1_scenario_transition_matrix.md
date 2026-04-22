# Phase 1 Scenario Transition Matrix

**Status**: Complete — gates Phase 3 implementation.  
**Date**: 2026-04-22  
**Source**: `public_api_symbol_inventory.md` (104 symbols), `public_api_full_surface_plan.md`, `issue_86_scenario_kernel.md`

---

## Canonical end-to-end user workflow

```mathematica
(* 1. Build a typed scenario from raw data *)
s = makeScenario[GetExampleData[12] /. {I1 -> 100, U1 -> 0}];

(* 2. Validate and complete *)
validateScenario[s]          (* returns True or issues a message *)
s = completeScenario[s]      (* fills defaults, normalizes *)

(* 3. Compile and solve *)
r = SolveMFG[s];

(* 4. Inspect *)
IsFeasible[r]
SolutionFlowPlot[r]
```

Backward-compatible path (existing scripts):

```mathematica
data = GetExampleData[12] /. {I1 -> 100, U1 -> 0};
d2e  = DataToEquations[data];
r    = CriticalCongestionSolver[d2e];
IsFeasible[r]
```

Both paths must work through at least Phase 6.

---

## Validation boundary policy

| Where | What is caught | What is NOT caught |
|---|---|---|
| `makeScenario` / `validateScenario` | Missing required keys, wrong value types, structural inconsistencies (e.g. adjacency matrix size mismatch), symbolic parameters not yet substituted | Numeric infeasibility, solver convergence |
| `DataToEquations` | Dimensionality / topology errors that produce malformed equation systems | Solver-level infeasibility |
| Solver (`CriticalCongestionSolver`, etc.) | Infeasible flow conditions, non-negative constraint violations | Nothing — solver returns structured result |

Rule: **fail fast at the scenario boundary** (makeScenario/validateScenario) for structural problems; let the solver report infeasibility as a typed result, not an exception.

---

## Quick-start outline for README rewrite

1. What is MFGraphs? (one paragraph)
2. Installation and loading
3. Minimal example: three lines (makeScenario → SolveMFG → IsFeasible)
4. Built-in examples: GetExampleData key list
5. Inspecting results: IsFeasible, SolutionFlowPlot, ExitFlowPlot
6. Advanced: direct solver access (CriticalCongestionSolver, DataToEquations)
7. Configuration globals ($MFGraphsVerbose, ClearSolveCache)
8. Links: API_REFERENCE.md, issue tracker

---

## Symbol classification

Decision codes:
- **core** — required for the scenario-first workflow; must have usage, tests, and docs
- **advanced-retained** — intentionally public but not in the basic workflow; must have usage
- **compatibility** — deprecated alias; keep through Phase 6 with deprecation message; retire in Phase 7+
- **retire** — remove from BeginPackage export list (dead symbol or accidental export)
- **examples-only** — move to MFGraphs\`Examples\` sub-context in Phase 6; not in main namespace
- **internal-blocked** — move to Private or FictitiousPlayBackend namespace; blocked by issues #72–#75

### Core (14 symbols)

| Symbol | Scenario workflow role | Phase 3 action |
|---|---|---|
| `SolveMFG` | Primary solver entrypoint | Accept `scenario[...]` as first arg |
| `DataToEquations` | Equation compiler | Accept `scenario[...]` as first arg |
| `CriticalCongestionSolver` | Primary solver (alpha=1) | No change (already works on compiled d2e) |
| `IsFeasible` | Solution inspection | No change |
| `GetExampleData` | Data access for examples | No change; wrap with makeScenario at call sites |
| `makeScenario` | Scenario construction | No change |
| `validateScenario` | Scenario validation | No change |
| `completeScenario` | Scenario normalization | No change |
| `scenarioQ` | Type predicate | No change |
| `scenario` | Typed head | No change |
| `ScenarioData` | Key-based accessor for scenario fields | No change |
| `NetworkGraphPlot` | Graph visualization | No change |
| `SolutionFlowPlot` | Flow visualization | No change |
| `ExitFlowPlot` | Exit flow visualization | No change |

### Advanced-retained (41 symbols)

These are intentionally public but outside the basic workflow. They serve research use,
solver developers, or downstream tooling. All must have `::usage` and at least one test by
Phase 7. No structural change in Phase 3.

**Solver internals (exposed for research/paper reproducibility)**:
| Symbol | Rationale |
|---|---|
| `IsCriticalSolution` | Allows callers to distinguish critical vs. non-critical solutions |
| `MFGPreprocessing` | Used in research scripts; inner solver stage |
| `MFGSystemSolver` | Inner solver stage; needed by NonLinear solver chain |
| `DNFSolveStep` | DNF step exposed for debugging / custom solver wiring |
| `TripleClean` | Fixed-point simplification; useful standalone |
| `TripleStep` | Single step of TripleClean; used in tests |
| `SystemToTriple` | System decomposition helper |
| `GetKirchhoffMatrix` | Matrix access for custom downstream tooling |
| `GetKirchhoffLinearSystem` | Same; explicit linear system form |
| `ConsistentSwitchingCosts` | Validate switching cost data before solve |
| `IsSwitchingCostConsistent` | Single-cost predicate variant |
| `DNFReduce` | Core algorithm; research users call directly |
| `SubstituteSolution` | Apply solution rules |
| `DeduplicateByComplexity` | DNF post-processing |
| `ReduceDisjuncts` | DNF post-processing |
| `ClearSolveCache` | Essential for multi-problem workflows |
| `CheckFlowFeasibility` | Feasibility check outside solver |
| `BuildSolverComparisonData` | Comparison tooling for paper |
| `EncodeFlowAssociation` | Encode flow data to vector |
| `DecodeFlowVector` | Decode vector back to association |
| `ComputeSignedEdgeFlowsFast` | Fast Kirchhoff residual support |
| `ComputeKirchhoffResidualFast` | Fast residual check |
| `NumberVectorQ` | Input guard |
| `LookupAssociationValue` | Utility |
| `SelectFlowAssociation` | Filter flow data |
| `BuildBoundaryMassData` | Boundary condition helper |
| `BuildUtilityReductionResidualData` | Residual helper |
| `FlowSplitting` | Flow structure helper |
| `FlowGathering` | Flow structure helper |
| `AltFlowOp` | Complementarity operator |
| `IneqSwitch` | Switching inequality helper |
| `AltSwitch` | Complementarity condition helper |
| `AssociationValue` | Graphics utility |
| `NetEdgeFlows` | Edge flow extraction |
| `NetworkVisualData` | Visual data builder |
| `FlowStyleDirective` | Style helper |

**Runtime configuration (must remain public)**:
| Symbol | Rationale |
|---|---|
| `alpha` | Congestion exponent — used by solver; user sets this for nonlinear cases |
| `Cost` | Congestion cost function — user-overridable |
| `j` | Flow variable head — appears in solutions |
| `u` | Utility variable head — appears in solutions |
| `z` | Vertex potential head — appears in solutions |
| `$MFGraphsVerbose` | Logging control |
| `$MFGraphsParallelThreshold` | Parallelism threshold |
| `$MFGraphsParallelReady` | Parallel kernel flag |
| `$MFGraphsParallelLaunchThreshold` | Launch threshold |
| `MFGPrint` | Verbose output helper |
| `MFGPrintTemporary` | Temporary output helper |
| `EnsureParallelKernels` | Kernel init |
| `MFGParallelMap` | Parallel dispatch |

**Note**: `alpha`, `j`, `u`, `z` are symbolic heads that appear in solution associations. They are "advanced-retained" not "core" because users only encounter them when inspecting raw solution payloads or building custom solvers — the scenario-first workflow insulates most users from them.

### Compatibility aliases (5 symbols)

Keep through Phase 6 with deprecation `::usage` text and `Message["Deprecated"]` on first call. Remove in Phase 7 or later.

| Symbol | Alias for | Action |
|---|---|---|
| `Data2Equations` | `DataToEquations` | Add deprecation message |
| `DataG` | `GetExampleData` | Add deprecation message |
| `FinalStep` | `DNFSolveStep` | Add deprecation message |
| `RemoveDuplicates` | `DeduplicateByComplexity` | Add deprecation message |
| `ReplaceSolution` | `SubstituteSolution` | Add deprecation message |

### Internal-blocked (9 symbols) — move to Private in Phase 6

These are Fictitious Play backend symbols blocked from public exposure by issues #72–#75.
In Phase 6, remove from BeginPackage; define with `Begin["`Private`"]` in `FictitiousPlayBackend.wl`.

| Symbol | Blocking issues |
|---|---|
| `SolveCriticalFictitiousPlayBackend` | #72, #73, #74, #75 |
| `SolveCriticalJFirstBackend` | #72, #73, #74, #75 |
| `SolveCriticalJFirstUtilities` | #72, #73, #74, #75 |
| `BuildFeasibleFlowSeed` | #72, #73 |
| `ExtractBellmanPotentials` | #72, #73 |
| `BuildSoftPolicyAndPropagate` | #73, #74 |
| `ClassifyAndCheckStability` | #73, #74 |
| `BuildCriticalQuadraticObjective` | #72, #75 |
| `UseQuadraticCriticalBackendQ` | #72, #75 |

### Orphaned monotone helpers (5 symbols) — retire in Phase 6

These are defined in `MFGraphs.wl` Private section but have no calling solver on
`codex/reduce-test-redundancy` (Monotone.wl was removed). They should move to
`rollback/pre-critical-only-surface` and be retired from master until the
Monotone solver is re-exposed as a `Method` option of `SolveMFG`.

| Symbol | Status |
|---|---|
| `BuildMonotoneStateData` | Orphaned; no caller on this branch |
| `BuildMonotoneValueSystem` | Orphaned; no caller on this branch |
| `BuildReducedKirchhoffCoordinates` | Orphaned; no caller on this branch |
| `BuildMonotonePairCostAssociation` | Orphaned; no caller on this branch |
| `MonotoneVariableFieldValue` | Orphaned; no caller on this branch |

**Phase 6 action**: remove `::usage` from BeginPackage; remove definitions from
`MFGraphs.wl` Private. Move definitions to a future `Monotone.wl` when it is
re-introduced.

### Examples-only parameters (22 symbols) — move namespace in Phase 6

These are symbolic parameter names exported alongside example data. They appear in
solution values when a user forgets to substitute, which can cause confusing output.

**Phase 6 action**: move `::usage` declarations and `BeginPackage` entries from the
root `MFGraphs\`` context to `MFGraphs\`Examples\``. Existing scripts using
`/. {I1 -> ..., U1 -> ...}` will continue to work if `GetExampleData` returns
symbols from the examples context.

| Symbol | Type |
|---|---|
| `I1`, `I2`, `I3` | Entry flow |
| `U1`, `U2`, `U3` | Terminal cost |
| `S1`–`S16` | Switching cost |

### Retire immediately (2 symbols)

These are accidental exports confirmed absent from the codebase by grep.

| Symbol | Reason |
|---|---|
| `alpha$` | Not defined anywhere; packaging artefact |
| `j$` | Not defined anywhere; packaging artefact |

**Phase 2 status**: both were already confirmed absent. Gate satisfied.

---

## Symbol count summary

| Decision | Count |
|---|---|
| core | 14 |
| advanced-retained | 41 |
| compatibility | 5 |
| internal-blocked | 9 |
| orphaned (retire in Phase 6) | 5 |
| examples-only (move namespace in Phase 6) | 22 |
| retire immediately | 2 |
| **Total** | **98** |

**Note**: The inventory listed 104 symbols but `alpha$` and `j$` were already absent,
and `$MFGraphsParallelLaunchThreshold` is newly added (not in original 104 count).
The effective live exported count is ~96. Exact count will be re-verified with
`Length[Names["MFGraphs\`*"]]` at Phase 6 gate.

---

## Phase 3 pre-requisites (what this matrix unlocks)

1. **`DataToEquations[s_scenario]`**: reads `s["Data"]`, passes Association to existing compiler.
2. **`SolveMFG[s_scenario]`**: calls `DataToEquations[s]` then dispatches to solver per `s["Method"]` or option.
3. Both must preserve backward-compat overloads `DataToEquations[assoc_Association]` and `SolveMFG[assoc_Association]`.
4. Validation failures from `validateScenario` must fire before solver entry — not deep inside symbolic reduction.

No symbol additions or removals are required for Phase 3. The scenario-first path is an additional overload, not a replacement.
