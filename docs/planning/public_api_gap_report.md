# MFGraphs Phase 1 public API gap report

This report summarizes the main gaps between the current exported symbol surface and the desired end state of a fully public, fully documented, and trusted MFGraphs API.

## Summary

| Metric | Count |
|---|---:|
| Exported symbols in current surface | 104 |
| Symbols with no direct test reference | 60 |
| Compatibility or deprecation symbols | 5 |
| Runtime globals exported as public symbols | 2 |
| Example-parameter symbols exported as public symbols | 22 |
| Suspicious or accidental exports | 2 |

## Highest-priority gaps

| Gap | Why it matters | Recommended Phase 2 action |
|---|---|---|
| Accidental exports such as `alpha$` and `j$` | These weaken confidence in the API boundary and pollute generated docs | Fix package export hygiene before broadening the public surface |
| Broad advanced-helper surface with sparse direct tests | Publicization without verification would reduce trust instead of increasing it | Add a symbol-to-test coverage matrix and targeted low-level unit tests |
| Runtime globals exposed without explicit stability policy | Users cannot tell whether these are supported controls or incidental internals | Decide which globals remain public and document the rest as internal controls |
| Example-parameter symbols exported beside functional API symbols | These clutter the generated reference and blur the distinction between examples and supported functions | Decide whether parameters stay exported, move to examples-only docs, or get grouped specially |
| Compatibility aliases mixed into the main public surface | This makes the API look larger and less coherent than it really is | Mark them consistently in docs and consider dedicated compatibility sections |

## Symbols with no direct test reference

| Symbol | Category | Decision |
|---|---|---|
| `$MFGraphsParallelReady` | runtime-global | review |
| `$MFGraphsVerbose` | runtime-global | review |
| `AltFlowOp` | advanced | review |
| `AltSwitch` | advanced | review |
| `AssociationValue` | advanced | review |
| `BuildBoundaryMassData` | advanced | review |
| `BuildCriticalQuadraticObjective` | advanced | review |
| `BuildFeasibleFlowSeed` | advanced | review |
| `BuildMonotonePairCostAssociation` | advanced | review |
| `BuildMonotoneStateData` | advanced | review |
| `BuildMonotoneValueSystem` | advanced | review |
| `BuildReducedKirchhoffCoordinates` | advanced | review |
| `BuildSoftPolicyAndPropagate` | advanced | review |
| `BuildSolverComparisonData` | advanced | review |
| `BuildUtilityReductionResidualData` | advanced | review |
| `CheckFlowFeasibility` | advanced | review |
| `ClassifyAndCheckStability` | advanced | review |
| `ClearSolveCache` | advanced | review |
| `ComputeKirchhoffResidualFast` | advanced | review |
| `ComputeSignedEdgeFlowsFast` | advanced | review |
| `ConsistentSwitchingCosts` | advanced | review |
| `DecodeFlowVector` | advanced | review |
| `DNFReduce` | advanced | review |
| `EncodeFlowAssociation` | advanced | review |
| `EnsureParallelKernels` | advanced | review |
| `ExtractBellmanPotentials` | advanced | review |
| `FlowGathering` | advanced | review |
| `FlowSplitting` | advanced | review |
| `FlowStyleDirective` | advanced | review |
| `IneqSwitch` | advanced | review |
| `IsSwitchingCostConsistent` | advanced | review |
| `LookupAssociationValue` | advanced | review |
| `MFGParallelMap` | advanced | review |
| `MFGPrint` | advanced | review |
| `MFGPrintTemporary` | advanced | review |
| `MonotoneVariableFieldValue` | advanced | review |
| `NetEdgeFlows` | advanced | review |
| `NetworkVisualData` | advanced | review |
| `NumberVectorQ` | advanced | review |
| `ReduceDisjuncts` | advanced | review |
| `SelectFlowAssociation` | advanced | review |
| `SolveCriticalJFirstBackend` | advanced | review |
| `SolveCriticalJFirstUtilities` | advanced | review |
| `SystemToTriple` | advanced | review |
| `TripleStep` | advanced | review |
| `UseQuadraticCriticalBackendQ` | advanced | review |
| `z` | symbolic-head | keep |

## Suspicious exports to fix first

| Symbol | Reason |
|---|---|
| `alpha$` | Exported in the generated API surface but not suitable as a stable public symbol |
| `j$` | Exported in the generated API surface but not suitable as a stable public symbol |

## Recommendation

Before making additional helpers intentionally public, Phase 2 should first normalize the exported surface: remove accidental exports, label compatibility symbols clearly, and define a stable classification rule for globals, symbolic heads, and example parameters. Only then should the repository expand usage contracts and per-symbol verification.
