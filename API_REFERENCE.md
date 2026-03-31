# MFGraphs API Reference

This document is automatically generated from the Wolfram Language package usage metadata.

## alpha

alpha[edge] is the edge-dependent congestion exponent used in the Hamiltonian. The default is 1.

## alphaFun

MFGraphs`alphaFun

## AltFlowOp

AltFlowOp[j][list] returns the alternative: j@@list ==0 || j@@Reverse@list ==0.

## AltSwitch

AltSwitch[j, u, switchingCosts][v, e1, e2] returns the complementarity condition:
(j[v, e1, e2] == 0) || (u[v, e1] == u[e2, e1] + switchingCosts[{v, e1, e2}])

## args

MFGraphs`args

## args$

MFGraphs`args$

## CachedGradientProjection

CachedGradientProjection[x, KM, dim, At, cache] is a version of GradientProjection
that caches the PseudoInverse result and reuses it when x has not changed significantly.
cache must be a symbol holding a 1-element list {Null} or {<|"x" -> ..., "pi" -> ...|>}.
The function has the HoldAll attribute so that cache is passed by reference.

## ClearSolveCache

ClearSolveCache[] clears the internal caches used by CachedSolve and CachedReduce.
Call this between different problem instances to prevent stale results.

## ConsistentSwitchingCosts

ConsistentSwitchingCosts[switchingcosts][{a,b,c}->S]
returns True if S, the cost of switching from the edge ab to cb, is smaller than any other combination,
such as, ab to bd and then from db to bc.
Returns the condition for this switching cost to satisfy the triangle inequality when S, and the other
switching costs too, does not have a numerical value.

## CriticalCongestionSolver

CriticalCongestionSolver[Eqs] returns Eqs with an association "AssoCritical" with rules to
the solution to the critical congestion case. Option: "ReturnShape" (default "Legacy"); use "Standard" to add normalized solver-result keys.

## Data2Equations

DataToEquations[Data] returns the equations, inequalities, and alternatives associated to the Data. 

## DataG

GetExampleData[n] returns an Association with keys "Vertices List", "Adjacency Matrix",
"Entrance Vertices and Flows", "Exit Vertices and Terminal Costs", and "Switching Costs"
from the test case test[n].
Supports test cases with 5 fields (basic), 6 fields (+cost function 'a'), or 7 fields (+alpha).

## DataToEquations

DataToEquations[Data] returns the equations, inequalities, and alternatives associated to the Data. 

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

## edge

MFGraphs`edge

## edge$

MFGraphs`edge$

## EnsureParallelKernels

EnsureParallelKernels[] launches parallel subkernels if none are running.

## expr

MFGraphs`expr

## f

MFGraphs`f

## FastIntegratedMass

FastIntegratedMass[interpM, j, edge] computes IntegratedMass using a precomputed interpolation of M.
interpM should be the result of PrecomputeM.

## FinalStep

DNFSolveStep[{EE,NN,OR}, rules] takes a grouped system and some Association of rules (a partial solution). It returns the result of applying DNFReduce.

## FlowGathering

FlowGathering[auxTriples_List][x_] returns the triples that end with x.

## FlowSplitting

FlowSplitting[AT][UndirectedEdge[a, b]] returns the splitting that start with {a,b}.

## g

g[m, edge] is the edge-dependent interaction term as a function of density m. The default is -1/m^2.

## GetExampleData

GetExampleData[n] returns an Association with keys "Vertices List", "Adjacency Matrix",
"Entrance Vertices and Flows", "Exit Vertices and Terminal Costs", and "Switching Costs"
from the test case test[n].
Supports test cases with 5 fields (basic), 6 fields (+cost function 'a'), or 7 fields (+alpha).

## GetKirchhoffLinearSystem

GetKirchhoffLinearSystem[d2e] returns the entry current vector, Kirchhoff matrix, and the variables in the order corresponding to the Kirchhoff matrix.

## GetKirchhoffMatrix

GetKirchhoffMatrix[d2e] returns the entry current vector, Kirchhoff matrix, (critical congestion) cost function placeholder, and the variables in the order corresponding to the Kirchhoff matrix. The third slot is retained for backward compatibility and should not be used by new code.

## GetTimeDependentExampleData

GetTimeDependentExampleData[key] returns an Association with all standard network keys
plus time-dependent keys ("Time Horizon", "Time Steps", etc.).
Available keys: "TD-1" through "TD-5".

## gFun

MFGraphs`gFun

## GradientProjection

GradientProjection[x, A, dim, At] returns the projected gradient operator matrix.
This is InverseHessian[x] . (I - At . PseudoInverse[A . InverseHessian[x] . At] . A . InverseHessian[x]).

## Hess

Hess[j] returns a numeric diagonal matrix with the reciprocals of the elements in the (numeric) vector j.

## HessianSandwich

HessianSandwich[j, A, At] returns the product A . InverseHessian[j] . At.

## I1

Symbolic parameter: input flow at entrance 1.

## I2

Symbolic parameter: input flow at entrance 2.

## I3

Symbolic parameter: input flow at entrance 3.

## IneqSwitch

IneqSwitch[u, switchingCosts][v, e1, e2] returns the optimality condition at the vertex v related to switching from e1 to e2. Namely,
u[v, e1] <= u[v, e2] + switchingCosts[{v, e1, e2}]

## InverseHessian

InverseHessian[j] returns a diagonal matrix with the (numeric) vector j. This is the inverse of the matrix Hess[j].

## IsFeasible

IsFeasible[result] returns True if a solver result is feasible, checking legacy "Status" or standardized "Feasibility" keys.

## IsNonLinearSolution

IsNonLinearSolution[Eqs] extracts AssoNonCritical and checks equations and inequalities. The right and left hand sides of the nonlinear equations are shown with the sup-norm of the difference.

## IsSwitchingCostConsistent

IsSwitchingCostConsistent[List of switching costs] is True if all switching costs satisfy the triangle inequality. If some switching costs are symbolic, then it returns the consistency conditions.

## IsTimeDependentQ

IsTimeDependentQ[data] returns True if the data Association contains time-dependent keys
(specifically, "Time Horizon").

## list

MFGraphs`list

## m

MFGraphs`m

## MFGParallelMap

MFGParallelMap[f, list] applies f to each element of list, using ParallelMap
when Length[list] >= $MFGraphsParallelThreshold, otherwise Map.

## MFGPreprocessing

MFGPreprocessing[Eqs] returns the association Eqs with the preliminary solution "InitRules" and corresponding 'reduced' "NewSystem".

## MFGPrint

MFGPrint[args___] prints args only when $MFGraphsVerbose is True.

## MFGPrintTemporary

MFGPrintTemporary[args___] prints a temporary message only when $MFGraphsVerbose is True.

## MFGSystemSolver

MFGSystemSolver[Eqs][edgeEquations] returns the
association with rules to the solution

## MonotoneSolver

MonotoneSolver[d2e] solves the MFG problem using the monotone operator method.
Options: "TimeSteps" (default 100), "UseCachedGradient" (default True), "ReturnShape" (default "Legacy"; use "Standard" for a normalized solver-result association), "PotentialFunction", "CongestionExponentFunction", and "InteractionFunction" (default Automatic, meaning use the current global MFGraphs` definitions of V, alpha, and g).

## MonotoneSolverFromData

MonotoneSolverFromData[Data] solves the MFG problem from raw Data using the monotone operator method.
Options: "TimeSteps" (default 100), "UseCachedGradient" (default True), "ReturnShape" (default "Legacy"; use "Standard" for a normalized solver-result association), "PotentialFunction", "CongestionExponentFunction", and "InteractionFunction" (default Automatic, meaning use the current global MFGraphs` definitions of V, alpha, and g).

## MonotoneSolverODE

MonotoneSolverODE[x0, KM, jj, cc] solves the gradient flow ODE from initial condition x0.
Options: "TimeSteps" (default 100), "UseCachedGradient" (default True).

## m$

MFGraphs`m$

## NonLinearSolver

NonLinearSolver[Eqs] takes an association resulting from DataToEquations and returns an approximation to the solution of the non-critical congestion case with alpha = value, specified by alpha[edge_] := value.
Options: "MaxIterations" (default 15), "Tolerance" (default 0), "ReturnShape" (default "Legacy"; use "Standard" to add normalized solver-result keys), "PotentialFunction", "CongestionExponentFunction", and "InteractionFunction" (default Automatic, meaning use the current global MFGraphs` definitions of V, alpha, and g). When Tolerance > 0, iteration stops early when the infinity-norm change in flow variables between consecutive steps falls below the given tolerance.

## NumberMatrixQ

NumberMatrixQ[A] returns True if the elements of the matrix A are numeric.

## NumberVectorQ

NumberVectorQ[j] returns True if the vector j is numeric.

## p

MFGraphs`p

## PrecomputeM

PrecomputeM[jMin, jMax, edge, nPoints] precomputes M[j,x,edge] on a grid and returns
an InterpolatingFunction. This avoids per-point FindRoot calls during NIntegrate.
nPoints controls the grid resolution (default 50).

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

## SubstituteSolution

SubstituteSolution[xp, sol] substitutes the Rule sol into the expression xp.
If the result has Head And, it simplifies the first conjunct.
Otherwise, it simplifies the whole expression.

## SystemToTriple

SystemToTriple[sys] returns a triple {equalities, inequalities, alternatives} from sys.
The input should be a system of equations, inequalities and (simple) alternatives.

## TimeDependentSolver

TimeDependentSolver[data] solves the time-dependent MFG problem on a network.
The data Association must include "Time Horizon" and optionally
"Time Steps", "Initial Mass Distribution", "Terminal Cost Function",
"Time Dependent Entrance Flows", and "Time Dependent Switching Costs".
Options: "MaxOuterIterations" (default 20), "Tolerance" (default 10^-6),
"SpatialSolverIterations" (default 0; when > 0, runs nonlinear iterations per step),
"ReturnShape" (default "Legacy"; use "Standard" for normalized output).

## TripleClean

TripleClean[{{EE,NN,OR},Rules}] composes TripleStep until it reaches a fixed point, that is, {{True,NewNN,NewOR},NewRules} such that replacement of NewRules in NewNN and NewOR do not produce equalities.

## TripleStep

TripleStep[{{EE,NN,OR},Rules}] returns {{NewEE, NewNN, NewOR}, NewRules}, where NewRules contain the solutions to all the equalities found in the system after replacing Rules in {EE,NN,OR}.

## U

U[x, edge, Eqs, sol] computes the value function at position x on the given edge.

## U1

Symbolic parameter: terminal cost at exit 1.

## U2

Symbolic parameter: terminal cost at exit 2.

## U3

Symbolic parameter: terminal cost at exit 3.

## V

V[x, edge] is the along-edge potential term in the Hamiltonian. Lower values make locations on an edge more attractive to agents. The default is 0.

## ValidateTimeDependentData

ValidateTimeDependentData[data] checks that the required time-dependent keys are present
and consistent. Returns True or a failure message string.

## vFun

MFGraphs`vFun

## WithHamiltonianFunctions

WithHamiltonianFunctions[vFun, alphaFun, gFun, expr] temporarily overrides the public MFGraphs Hamiltonian ingredients V, alpha, and g while evaluating expr. Any argument may be Automatic to keep the current definition.

## x

MFGraphs`x

## xi

MFGraphs`xi

## x$

MFGraphs`x$

## $MFGraphsParallelReady

False is the symbol for the Boolean value false. 

## $MFGraphsVerbose

False is the symbol for the Boolean value false. 
