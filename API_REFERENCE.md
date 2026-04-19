# MFGraphs API Reference

This document is automatically generated from the Wolfram Language package usage metadata.

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

## NumberMatrixQ

NumberMatrixQ[A] returns True if the elements of the matrix A are numeric.

## NumberVectorQ

NumberVectorQ[j] returns True if the vector j is numeric.

## p

MFGraphs`p

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

## x

MFGraphs`x

## xi

MFGraphs`xi

## x$

MFGraphs`x$

## $MFGraphsParallelReady

False is the symbol for the Boolean value false. 

## $MFGraphsVerbose

True is the symbol for the Boolean value true. 
