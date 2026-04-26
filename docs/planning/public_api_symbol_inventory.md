# MFGraphs Phase 1 public API symbol inventory

This document inventories the currently exported symbols in the **MFGraphs public context** and classifies each symbol as part of the Phase 1 publicization review.

The current exported surface contains **104** symbols.

## Inventory table

| Symbol | Declared in | Category | Decision | Stability | Test hits | Notes |
|---|---|---|---|---|---:|---|
| `$MFGraphsParallelReady` | `MFGraphs/MFGraphs.wl` | runtime-global | review | review-needed | 0 | no direct test reference found; runtime global |
| `$MFGraphsVerbose` | `MFGraphs/MFGraphs.wl` | runtime-global | review | review-needed | 0 | no direct test reference found; runtime global |
| `alpha` | `MFGraphs/MFGraphs.wl` | symbolic-head | keep | stable-candidate | 1 | — |
| `alpha$` | — | accidental | fix-export | broken-export | 0 | appears accidental; no direct test reference found |
| `AltFlowOp` | `MFGraphs/DataToEquations.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `AltSwitch` | `MFGraphs/DataToEquations.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `AssociationValue` | `MFGraphs/Graphics.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `BuildBoundaryMassData` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `BuildCriticalQuadraticObjective` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `BuildFeasibleFlowSeed` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `BuildMonotonePairCostAssociation` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `BuildMonotoneStateData` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `BuildMonotoneValueSystem` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `BuildReducedKirchhoffCoordinates` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `BuildSoftPolicyAndPropagate` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `BuildSolverComparisonData` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `BuildUtilityReductionResidualData` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `CheckFlowFeasibility` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `ClassifyAndCheckStability` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `ClearSolveCache` | `MFGraphs/DNFReduce.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `completeScenario` | `MFGraphs/Scenario.wl` | core | keep | stable-candidate | 1 | — |
| `ComputeKirchhoffResidualFast` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `ComputeSignedEdgeFlowsFast` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `ConsistentSwitchingCosts` | `MFGraphs/DataToEquations.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `Cost` | `MFGraphs/MFGraphs.wl` | symbolic-head | keep | stable-candidate | 1 | — |
| `CriticalCongestionSolver` | `MFGraphs/MFGraphs.wl`, `MFGraphs/Solvers.wl` | core | keep | stable-candidate | 14 | duplicate usage declaration |
| `Data2Equations` | `MFGraphs/DataToEquations.wl` | compatibility | keep-with-deprecation | compatibility | 1 | — |
| `DataG` | `MFGraphs/Examples/ExamplesData.wl` | compatibility | review | compatibility | 1 | — |
| `DataToEquations` | `MFGraphs/DataToEquations.wl`, `MFGraphs/MFGraphs.wl` | core | keep | stable-candidate | 15 | duplicate usage declaration |
| `DecodeFlowVector` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `DeduplicateByComplexity` | `MFGraphs/DNFReduce.wl` | advanced | review | review-needed | 1 | — |
| `DNFReduce` | `MFGraphs/DNFReduce.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `DNFSolveStep` | `MFGraphs/DataToEquations.wl` | advanced | review | review-needed | 1 | — |
| `EncodeFlowAssociation` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `EnsureParallelKernels` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `ExitFlowPlot` | `MFGraphs/Graphics.wl` | core | keep | stable-candidate | 1 | — |
| `ExtractBellmanPotentials` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `FinalStep` | `MFGraphs/DataToEquations.wl` | compatibility | keep-with-deprecation | compatibility | 1 | — |
| `FlowGathering` | `MFGraphs/DataToEquations.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `FlowSplitting` | `MFGraphs/DataToEquations.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `FlowStyleDirective` | `MFGraphs/Graphics.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `GetExampleData` | `MFGraphs/Examples/ExamplesData.wl` | advanced | review | review-needed | 15 | — |
| `GetKirchhoffLinearSystem` | `MFGraphs/DataToEquations.wl` | advanced | review | review-needed | 1 | — |
| `GetKirchhoffMatrix` | `MFGraphs/DataToEquations.wl` | advanced | review | review-needed | 1 | — |
| `I1` | `MFGraphs/Examples/ExamplesData.wl` | example-parameter | review | review-needed | 15 | example parameter symbol |
| `I2` | `MFGraphs/Examples/ExamplesData.wl` | example-parameter | review | review-needed | 6 | example parameter symbol |
| `I3` | `MFGraphs/Examples/ExamplesData.wl` | example-parameter | review | review-needed | 0 | no direct test reference found; example parameter symbol |
| `IneqSwitch` | `MFGraphs/DataToEquations.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `IsCriticalSolution` | `MFGraphs/MFGraphs.wl`, `MFGraphs/Solvers.wl` | core | keep | stable-candidate | 3 | duplicate usage declaration |
| `IsFeasible` | `MFGraphs/MFGraphs.wl` | core | keep | stable-candidate | 8 | — |
| `IsSwitchingCostConsistent` | `MFGraphs/DataToEquations.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `j` | `MFGraphs/MFGraphs.wl` | symbolic-head | keep | stable-candidate | 6 | — |
| `j$` | — | accidental | fix-export | broken-export | 0 | appears accidental; no direct test reference found |
| `LookupAssociationValue` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `makeScenario` | `MFGraphs/Scenario.wl` | core | keep | stable-candidate | 1 | — |
| `MFGParallelMap` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `MFGPreprocessing` | `MFGraphs/DataToEquations.wl`, `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 1 | duplicate usage declaration |
| `MFGPrint` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `MFGPrintTemporary` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `MFGSystemSolver` | `MFGraphs/DataToEquations.wl`, `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 3 | duplicate usage declaration |
| `MonotoneVariableFieldValue` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `NetEdgeFlows` | `MFGraphs/Graphics.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `NetworkGraphPlot` | `MFGraphs/Graphics.wl` | core | keep | stable-candidate | 1 | — |
| `NetworkVisualData` | `MFGraphs/Graphics.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `NumberVectorQ` | `MFGraphs/DataToEquations.wl`, `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | duplicate usage declaration; no direct test reference found |
| `ReduceDisjuncts` | `MFGraphs/DNFReduce.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `RemoveDuplicates` | `MFGraphs/DNFReduce.wl` | compatibility | keep-with-deprecation | compatibility | 1 | — |
| `ReplaceSolution` | `MFGraphs/DNFReduce.wl` | compatibility | keep-with-deprecation | compatibility | 1 | — |
| `S1` | `MFGraphs/Examples/ExamplesData.wl` | example-parameter | review | review-needed | 4 | example parameter symbol |
| `S10` | `MFGraphs/Examples/ExamplesData.wl` | example-parameter | review | review-needed | 0 | no direct test reference found; example parameter symbol |
| `S11` | `MFGraphs/Examples/ExamplesData.wl` | example-parameter | review | review-needed | 0 | no direct test reference found; example parameter symbol |
| `S12` | `MFGraphs/Examples/ExamplesData.wl` | example-parameter | review | review-needed | 0 | no direct test reference found; example parameter symbol |
| `S13` | `MFGraphs/Examples/ExamplesData.wl` | example-parameter | review | review-needed | 0 | no direct test reference found; example parameter symbol |
| `S14` | `MFGraphs/Examples/ExamplesData.wl` | example-parameter | review | review-needed | 0 | no direct test reference found; example parameter symbol |
| `S15` | `MFGraphs/Examples/ExamplesData.wl` | example-parameter | review | review-needed | 0 | no direct test reference found; example parameter symbol |
| `S16` | `MFGraphs/Examples/ExamplesData.wl` | example-parameter | review | review-needed | 0 | no direct test reference found; example parameter symbol |
| `S2` | `MFGraphs/Examples/ExamplesData.wl` | example-parameter | review | review-needed | 3 | example parameter symbol |
| `S3` | `MFGraphs/Examples/ExamplesData.wl` | example-parameter | review | review-needed | 1 | example parameter symbol |
| `S4` | `MFGraphs/Examples/ExamplesData.wl` | example-parameter | review | review-needed | 1 | example parameter symbol |
| `S5` | `MFGraphs/Examples/ExamplesData.wl` | example-parameter | review | review-needed | 1 | example parameter symbol |
| `S6` | `MFGraphs/Examples/ExamplesData.wl` | example-parameter | review | review-needed | 1 | example parameter symbol |
| `S7` | `MFGraphs/Examples/ExamplesData.wl` | example-parameter | review | review-needed | 0 | no direct test reference found; example parameter symbol |
| `S8` | `MFGraphs/Examples/ExamplesData.wl` | example-parameter | review | review-needed | 0 | no direct test reference found; example parameter symbol |
| `S9` | `MFGraphs/Examples/ExamplesData.wl` | example-parameter | review | review-needed | 0 | no direct test reference found; example parameter symbol |
| `scenario` | `MFGraphs/Scenario.wl` | advanced | review | review-needed | 1 | — |
| `ScenarioData` | `MFGraphs/Scenario.wl` | advanced | keep | stable-candidate | 1 | — |
| `scenarioQ` | `MFGraphs/Scenario.wl` | advanced | keep | stable-candidate | 1 | — |
| `SelectFlowAssociation` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `SolutionFlowPlot` | `MFGraphs/Graphics.wl` | core | keep | stable-candidate | 1 | — |
| `SolveCriticalFictitiousPlayBackend` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 1 | — |
| `SolveCriticalJFirstBackend` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `SolveCriticalJFirstUtilities` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `SolveMFG` | `MFGraphs/MFGraphs.wl`, `MFGraphs/SolveMFGDispatch.wl` | core | keep | stable-candidate | 2 | duplicate usage declaration |
| `SubstituteSolution` | `MFGraphs/DNFReduce.wl` | advanced | review | review-needed | 1 | — |
| `SystemToTriple` | `MFGraphs/DataToEquations.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `TripleClean` | `MFGraphs/DataToEquations.wl` | advanced | review | review-needed | 1 | — |
| `TripleStep` | `MFGraphs/DataToEquations.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `u` | `MFGraphs/MFGraphs.wl` | symbolic-head | keep | stable-candidate | 2 | — |
| `U1` | `MFGraphs/Examples/ExamplesData.wl` | example-parameter | review | review-needed | 15 | example parameter symbol |
| `U2` | `MFGraphs/Examples/ExamplesData.wl` | example-parameter | review | review-needed | 11 | example parameter symbol |
| `U3` | `MFGraphs/Examples/ExamplesData.wl` | example-parameter | review | review-needed | 5 | example parameter symbol |
| `UseQuadraticCriticalBackendQ` | `MFGraphs/MFGraphs.wl` | advanced | review | review-needed | 0 | no direct test reference found |
| `validateScenario` | `MFGraphs/Scenario.wl` | core | keep | stable-candidate | 1 | — |
| `z` | `MFGraphs/MFGraphs.wl` | symbolic-head | keep | stable-candidate | 0 | no direct test reference found |

## Notable private implementation families

The list below is intentionally approximate. It highlights representative helper definitions visible after entering the private section of each source module, which makes them useful candidates for later wrapper extraction or publicization review.

| Module | Representative private helper names |
|---|---|
| `MFGraphs/DNFReduce.wl` | `CachedSolve`, `CachedReduce`, `FlattenConjuncts`, `SubsumptionPrune` |
| `MFGraphs/DataToEquations.wl` | `TransitionsAt`, `ComputeDistanceToExitAssociation`, `BuildSignedEdgeMatrixFromKirchhoff`, `ExtractDirectedEdgePairsFromAdjacencyMatrix`, `BuildCriticalDecouplingPartition`, `BuildNumericState`, `RoundValues`, `UtilityConjunctionTerms`, `CanonicalUtilityOrderKey`, `BuildUtilityReductionData`, `BuildPrunedSystem`, `BuildOraclePrunedSystem` |
| `MFGraphs/Graphics.wl` | `vertexPropertyMap` |
| `MFGraphs/MFGraphs.wl` | `Hess`, `InverseHessian`, `BuildCriticalQuadraticEdgeModel`, `MakeSolverResult`, `SolveMFGCompiledInputQ` |
| `MFGraphs/Solvers.wl` | `GraphDistanceHeuristicSafeQ`, `ResolveOneOr`, `ResolveOrByGraphDistance`, `DirectCriticalSolver`, `BuildUEquivalenceReduction`, `CriticalNumericBackendRequestedQ`, `CriticalNumericBackendEligibleQ`, `BuildCriticalLinearConstraints`, `BuildLinearSystemFromEqualities`, `LinearSolveCandidate`, `CriticalJFirstFailure`, `CriticalJFirstFailureReason` |

## Phase 1 observations

The exported surface is broader than the currently documented user narrative. It mixes core workflow functions, advanced solver helpers, runtime globals, symbolic heads, compatibility aliases, and several example-parameter symbols. It also contains suspicious exports such as `alpha$` and `j$`, which appear to be packaging artefacts rather than intentional API symbols.
