# MFGraphs API Reference

This document is automatically generated from the Wolfram Language package usage metadata.

## alpha

alpha[edge] is the congestion exponent for an edge. Default is 1.

## alpha$

MFGraphs`alpha$

## AltFlowOp

AltFlowOp[j][list] returns the alternative: j@@list ==0 || j@@Reverse@list ==0.

## AltSwitch

AltSwitch[j, u, switchingCosts][v, e1, e2] returns the complementarity condition:
(j[v, e1, e2] == 0) || (u[v, e1] == u[e2, e1] + switchingCosts[{v, e1, e2}])

## AssociationValue

AssociationValue[assoc, key, default] returns assoc[key] when key exists, otherwise default (Missing["NotAvailable"] by default).

## BuildBoundaryMassData

BuildBoundaryMassData[Eqs, flowAssoc] returns mass balance metrics.

## BuildCriticalQuadraticObjective

BuildCriticalQuadraticObjective[d2e] returns a quadratic approximation of the MFG objective.

## BuildFeasibleFlowSeed

BuildFeasibleFlowSeed[backendState] returns an initial feasible flow for FP.

## BuildMonotonePairCostAssociation

BuildMonotonePairCostAssociation[halfPairs, edgeList, q] returns signed edge costs.

## BuildMonotoneStateData

BuildMonotoneStateData[d2e] returns shared linear state for monotone-like solvers.

## BuildMonotoneValueSystem

BuildMonotoneValueSystem[d2e] returns a function to solve for node potentials.

## BuildReducedKirchhoffCoordinates

BuildReducedKirchhoffCoordinates[d2e, basePoint] returns reduced affine coordinates on the Kirchhoff manifold.

## BuildSoftPolicyAndPropagate

BuildSoftPolicyAndPropagate[backendState, flowState, potentialState] propagates policy.

## BuildSolverComparisonData

BuildSolverComparisonData[Eqs, solution] returns an association with comparison metrics like Kirchhoff residual and boundary mass balance.

## BuildUtilityReductionResidualData

BuildUtilityReductionResidualData[Eqs, solution] returns utility reduction metrics.

## CheckFlowFeasibility

CheckFlowFeasibility[assoc] returns "Feasible" if all flow variables in assoc are non-negative, and "Infeasible" otherwise.

## ClassifyAndCheckStability

ClassifyAndCheckStability[backendState, flowState, potentialState, history] checks stability.

## ClearSolveCache

ClearSolveCache[] clears the internal caches used by CachedSolve and CachedReduce.
Call this between different problem instances to prevent stale results.

## completeScenario

completeScenario[s] fills in derived fields: sets "contentHash" in the Identity block (SHA256 of the canonical Model string), and supplies default Benchmark values ("Tier" -> "core", "Timeout" -> 300) when missing. Returns a new scenario object.

## ComputeKirchhoffResidualFast

ComputeKirchhoffResidualFast[ns, vec] returns the maximum Kirchhoff residual.

## ComputeSignedEdgeFlowsFast

ComputeSignedEdgeFlowsFast[ns, vec] returns a vector of signed edge flows.

## ConsistentSwitchingCosts

ConsistentSwitchingCosts[switchingcosts][{a,b,c}->S]
returns True if S, the cost of switching from the edge ab to cb, is smaller than any other combination,
such as, ab to bd and then from db to bc.
Returns the condition for this switching cost to satisfy the triangle inequality when S, and the other
switching costs too, does not have a numerical value.

## Cost

Cost[m, edge] is the congestion cost function.

## CriticalCongestionSolver

CriticalCongestionSolver[Eqs] returns Eqs with an association "AssoCritical" with rules to
the solution to the critical congestion case. It returns a standardized association
containing solver metadata, feasibility, comparison fields, and the solver-specific
payload key "AssoCritical". Options: "ValidateSolution" (default True), "ValidationTolerance" (default $CriticalSolverTolerance), "ValidationVerbose" (default False), "SymbolicTimeLimit" (default 120., time budget in seconds for the symbolic pipeline), and "ExactMode" (default False; when True, skips numeric/direct/oracle paths and returns "SymbolicRegion" for undetermined variables).

## Data2Equations

Data2Equations[Data] is a deprecated compatibility wrapper for DataToEquations[Data]. It will be removed in the next release.

## DataG

GetExampleData[n] returns an Association with keys "Vertices List", "Adjacency Matrix",
"Entrance Vertices and Flows", "Exit Vertices and Terminal Costs", and "Switching Costs"
from the test case test[n].
Supports test cases with 5 fields (basic), 6 fields (+cost function 'a'), or 7 fields (+alpha).

## DataToEquations

DataToEquations[Data] returns the equations, inequalities, and alternatives associated to the Data. 

## DecodeFlowVector

DecodeFlowVector[ns, vec] returns a flow association from a packed numeric vector.

## DeduplicateByComplexity

DeduplicateByComplexity[xp] sorts (by SimplifyCount) and then DeleteDuplicates.
For example, (A && B) || (B && A) becomes (A && B) only after sorting.

## DNFReduce

DNFReduce[xp, sys] converts the system xp && sys into disjunctive normal form
by solving equalities, substituting, and reducing branches.
DNFReduce[xp, sys, elem] is the 3-argument form that processes a single element
from a conjunction: if elem is an equality, solve and substitute into xp and sys;
otherwise, absorb elem into xp and continue with sys.

## DNFSolveStep

DNFSolveStep[{EE,NN,OR}, rules] takes a grouped system and some Association of rules (a partial solution). It returns the result of applying DNFReduce.

## EncodeFlowAssociation

EncodeFlowAssociation[ns, assoc] returns a packed numeric vector from a flow association.

## EnsureParallelKernels

EnsureParallelKernels[] launches parallel subkernels if none are running.

## ExitFlowPlot

ExitFlowPlot[exitFlows] produces a bar chart of exit-flow totals from an association mapping exit vertices to numeric flow values. An optional second argument sets the plot title.

## ExtractBellmanPotentials

ExtractBellmanPotentials[backendState, flowState] extracts node potentials.

## FinalStep

DNFSolveStep[{EE,NN,OR}, rules] takes a grouped system and some Association of rules (a partial solution). It returns the result of applying DNFReduce.

## FlowGathering

FlowGathering[auxTriples_List][x_] returns the triples that end with x.

## FlowSplitting

FlowSplitting[AT][UndirectedEdge[a, b]] returns the splitting that start with {a,b}.

## FlowStyleDirective

FlowStyleDirective[flow, maxFlow] returns an edge style directive (color/thickness/opacity) scaled by flow magnitude and sign.

## GetExampleData

GetExampleData[n] returns an Association with keys "Vertices List", "Adjacency Matrix",
"Entrance Vertices and Flows", "Exit Vertices and Terminal Costs", and "Switching Costs"
from the test case test[n].
Supports test cases with 5 fields (basic), 6 fields (+cost function 'a'), or 7 fields (+alpha).

## GetKirchhoffLinearSystem

GetKirchhoffLinearSystem[d2e] returns the entry current vector, Kirchhoff matrix, and the variables in the order corresponding to the Kirchhoff matrix.

## GetKirchhoffMatrix

GetKirchhoffMatrix[d2e] returns the entry current vector, Kirchhoff matrix, (critical congestion) cost function placeholder, and the variables in the order corresponding to the Kirchhoff matrix. The third slot is retained for backward compatibility and should not be used by new code.

## I1

Symbolic parameter: input flow at entrance 1.

## I2

Symbolic parameter: input flow at entrance 2.

## I3

Symbolic parameter: input flow at entrance 3.

## IneqSwitch

IneqSwitch[u, switchingCosts][v, e1, e2] returns the optimality condition at the vertex v related to switching from e1 to e2. Namely,
u[v, e1] <= u[v, e2] + switchingCosts[{v, e1, e2}]

## IsCriticalSolution

IsCriticalSolution[Eqs] validates whether the critical-congestion solution stored in "AssoCritical" satisfies the full symbolic MFG constraint set (equalities, inequalities, alternatives/complementarity) and the critical EqGeneral residual. By default it returns True/False. Options: "Tolerance" (default $CriticalSolverTolerance), "Verbose" (default False), and "ReturnReport" (default False).

## IsFeasible

IsFeasible[result] returns True if the solver result indicates a feasible solution (all flows non-negative and constraints satisfied within tolerance).

## IsSwitchingCostConsistent

IsSwitchingCostConsistent[List of switching costs] is True if all switching costs satisfy the triangle inequality. If some switching costs are symbolic, then it returns the consistency conditions.

## j

j[v, e] or j[v, e1, e2] represents a flow variable.

## j$

MFGraphs`j$

## LookupAssociationValue

LookupAssociationValue[assoc, key, default] is a robust Lookup helper.

## makeScenario

makeScenario[assoc] constructs a typed scenario from a raw association. The input must contain a "Model" key whose value is a network topology association accepted by DataToEquations (required keys: "Vertices List", "Adjacency Matrix", "Entrance Vertices and Flows", "Exit Vertices and Terminal Costs", "Switching Costs"). Optional keys: "Data" (parameter substitution rules), "Identity" (name, version), "Benchmark" (tier, timeout), "Visualization", "Inheritance". Returns a scenario[...] object on success or Failure[...] on error.

## MFGParallelMap

MFGParallelMap[f, list] applies f to each element of list, using ParallelMap
when Length[list] >= $MFGraphsParallelThreshold (and kernels are already running)
or Length[list] >= $MFGraphsParallelLaunchThreshold (when kernels must be launched).

## MFGPreprocessing

MFGPreprocessing[Eqs] returns the association Eqs with the preliminary solution "InitRules" and corresponding 'reduced' "NewSystem".

## MFGPrint

MFGPrint[args___] prints args only when $MFGraphsVerbose is True.

## MFGPrintTemporary

MFGPrintTemporary[args___] prints a temporary message only when $MFGraphsVerbose is True.

## MFGSystemSolver

MFGSystemSolver[Eqs][edgeEquations] returns an Association <|"Solution" -> assoc_or_Null, "UnresolvedConstraints" -> expr_or_None|>. "Solution" is Null on failure; "UnresolvedConstraints" is non-None when flows were determined but u-variables remain symbolically underdetermined.

## MonotoneVariableFieldValue

MonotoneVariableFieldValue[var, values, switching] returns the field value for a variable.

## NetEdgeFlows

NetEdgeFlows[d2e, solution, pairs] returns the net signed flow on each requested edge pair after applying balance and boundary flow rules. pairs defaults to all model edges.

## NetworkGraphPlot

NetworkGraphPlot[d2e] plots the network structure contained in the standardized DataToEquations association d2e. An optional second argument sets the plot title. The function is intended for workbook and package-level visualization of the directed network topology.

## NetworkVisualData

NetworkVisualData[d2e] builds reusable graph layout and vertex styling metadata used by plotting helpers.

## NumberVectorQ

NumberVectorQ[j] returns True if the vector j is numeric.

## ReduceDisjuncts

ReduceDisjuncts[expr] reduces the number of disjuncts in a DNF expression by
removing subsumed branches. For small expressions (<= $ReduceDisjunctsThreshold
disjuncts) uses Reduce directly; for larger expressions uses subsumption pruning
followed by Reduce on the survivors.

## RemoveDuplicates

DeduplicateByComplexity[xp] sorts (by SimplifyCount) and then DeleteDuplicates.
For example, (A && B) || (B && A) becomes (A && B) only after sorting.

## ReplaceSolution

SubstituteSolution[xp, sol] substitutes the Rule sol into the expression xp.
If the result has Head And, it simplifies the first conjunct.
Otherwise, it simplifies the whole expression.

## S1

Symbolic parameter: switching cost 1.

## S10

Symbolic parameter: switching cost 10.

## S11

Symbolic parameter: switching cost 11.

## S12

Symbolic parameter: switching cost 12.

## S13

Symbolic parameter: switching cost 13.

## S14

Symbolic parameter: switching cost 14.

## S15

Symbolic parameter: switching cost 15.

## S16

Symbolic parameter: switching cost 16.

## S2

Symbolic parameter: switching cost 2.

## S3

Symbolic parameter: switching cost 3.

## S4

Symbolic parameter: switching cost 4.

## S5

Symbolic parameter: switching cost 5.

## S6

Symbolic parameter: switching cost 6.

## S7

Symbolic parameter: switching cost 7.

## S8

Symbolic parameter: switching cost 8.

## S9

Symbolic parameter: switching cost 9.

## scenario

scenario[assoc] is the typed head for a MFGraphs scenario object. Use makeScenario to construct one; use ScenarioData to access keys.

## ScenarioData

ScenarioData[s, key] returns the value associated with key in the scenario s, or Missing["KeyAbsent", key] if absent. ScenarioData[s] returns the underlying Association.

## scenarioQ

scenarioQ[x] returns True if x is a typed scenario[assoc_Association] object, False otherwise.

## SelectFlowAssociation

SelectFlowAssociation[assoc] returns a sub-association with only flow (j) keys.

## SolutionFlowPlot

SolutionFlowPlot[d2e, solution] plots the network with edge labels and edge styling determined by the net edge flows implied by solution. An optional third argument sets the plot title. The input d2e should be the standardized DataToEquations association and solution should be an association of replacement rules or solved values compatible with its symbolic flow expressions.

## SolveCriticalFictitiousPlayBackend

SolveCriticalFictitiousPlayBackend[backendState] runs the FP iterative solver.

## SolveCriticalJFirstBackend

SolveCriticalJFirstBackend[Eqs, opts] runs the j-first numeric strategy.

## SolveCriticalJFirstUtilities

SolveCriticalJFirstUtilities[eqs, flowAssoc, uVars, tol] recovers node potentials from flows.

## SolveMFG

SolveMFG[input, opts] solves the critical congestion MFG and returns its standardized result association.

input can be either raw data (association accepted by DataToEquations) or an already compiled equation association.

## SubstituteSolution

SubstituteSolution[xp, sol] substitutes the Rule sol into the expression xp.
If the result has Head And, it simplifies the first conjunct.
Otherwise, it simplifies the whole expression.

## SystemToTriple

SystemToTriple[sys] returns a triple {equalities, inequalities, alternatives} from sys.
The input should be a system of equations, inequalities and (simple) alternatives.

## TripleClean

TripleClean[{{EE,NN,OR},Rules}] composes TripleStep until it reaches a fixed point, that is, {{True,NewNN,NewOR},NewRules} such that replacement of NewRules in NewNN and NewOR do not produce equalities.

## TripleStep

TripleStep[{{EE,NN,OR},Rules}] returns {{NewEE, NewNN, NewOR}, NewRules}, where NewRules contain the solutions to all the equalities found in the system after replacing Rules in {EE,NN,OR}.

## u

u[v, e] represents a utility/potential variable.

## U1

Symbolic parameter: terminal cost at exit 1.

## U2

Symbolic parameter: terminal cost at exit 2.

## U3

Symbolic parameter: terminal cost at exit 3.

## UseQuadraticCriticalBackendQ

UseQuadraticCriticalBackendQ[d2e] returns True if the case is eligible for a quadratic shortcut.

## validateScenario

validateScenario[s] checks that the scenario s has all required Model keys and that the Model value is an Association. Returns s unchanged on success, or Failure["ScenarioValidation", <|"Message" -> msg, "MissingKeys" -> {...}|>] on failure.

## z

z[v] represents a vertex potential variable.

## $MFGraphsParallelReady

False is the symbol for the Boolean value false. 

## $MFGraphsVerbose

False is the symbol for the Boolean value false. 
