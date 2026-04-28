# MFGraphs API Reference

This document is automatically generated from the Wolfram Language package usage metadata.

## alpha

alpha[edge] is the congestion exponent for an edge. Default is 1.

## Cost

Cost[m, edge] is the congestion cost function.

## ensureParallelKernels

ensureParallelKernels[] launches parallel subkernels if none are running.

## j

j[v, e] or j[v, e1, e2] represents a flow variable.

## mfgPrint

mfgPrint[args___] prints args only when $MFGraphsVerbose is True.

## mfgPrintTemporary

mfgPrintTemporary[args___] prints a temporary message only when $MFGraphsVerbose is True.

## u

u[v, e] represents a value-function variable. u[r, i] is the value at vertex v_i on edge e_{r,i}.

## z

z[v] represents a vertex potential variable.

## buildAuxiliaryTopology

buildAuxiliaryTopology[model] returns an association with the auxiliary graph and metadata derived from the raw model.

## buildAuxTriples

buildAuxTriples[auxGraph] returns the list of all possible {v_in, v_mid, v_out} transitions (triples) in the graph.

## completeScenario

completeScenario[s] fills in derived fields and supplies default Benchmark values ("Tier" -> "core", "Timeout" -> 300) when missing. If cached topology is absent, it warns and rebuilds topology from the Model before completing. Returns a new scenario object.

## deriveAuxPairs

deriveAuxPairs[topology] returns the list of all directed edge pairs {u, v} in the auxiliary graph (including entry/exit and reversed graph edges).

## makeScenario

makeScenario[assoc] constructs a typed scenario from a raw association. The input must contain a "Model" key whose value is a network topology association (required keys: "Vertices", "Adjacency", "Entries", "Exits", "Switching"). If "Model" contains a Wolfram Graph under key "Graph", missing "Vertices" and/or "Adjacency" are derived automatically. Optional keys: "Hamiltonian" (<|"Alpha" -> a, "V" -> v, "G" -> g, "EdgeAlpha" -> <|{u,v} -> a_uv, ...|>, "EdgeV" -> <|{u,v} -> v_uv, ...|>, "EdgeG" -> <|{u,v} -> g_uv, ...|>|>), "Identity" (name, version), "Benchmark" (tier, timeout), "Visualization", "Inheritance". Default Hamiltonian is Alpha=1 and V=0 on all edges, with G[z]=-1/z (overridable globally and per edge). Boundary values must be numeric; switching-cost values must be numeric or Infinity. Returns a scenario[...] object on success or Failure[...] on error.

## scenario

scenario[assoc] is the typed head for a MFGraphs scenario object. Use makeScenario to construct one; use scenarioData to access keys.

## scenarioData

scenarioData[s, key] returns the value associated with key in the scenario s, or Missing["KeyAbsent", key] if absent. scenarioData[s] returns the underlying Association.

## scenarioQ

scenarioQ[x] returns True if x is a typed scenario[assoc_Association] object, False otherwise.

## validateScenario

validateScenario[s] checks that the scenario s has all required Model keys and that the Model value is an Association. Returns s unchanged on success, or Failure["ScenarioValidation", <|"Message" -> msg, "MissingKeys" -> {...}|>] on failure.

## amScenario

amScenario[vl, am, entries, exits] creates a scenario from an explicit vertices list vl and adjacency matrix am. Vertex labels in vl must be positive integers; non-integer labels are not accepted. Optional: sc, alpha, V, g.

## cycleScenario

cycleScenario[n, entries, exits] creates a scenario on a directed n-cycle (1->2->...->n->1), vertices 1..n. Optional: sc, alpha, V, g.

## getExampleScenario

getExampleScenario[n] returns a 6-arg factory Function[{entries,exits,sc,alpha,V,g}, scenario[...]] for built-in example n. Topology is baked in; all parameters are caller-supplied. getExampleScenario[n, entries, exits] calls the factory using the canonical switching costs for that case (sc=Automatic resolves via $CaseDefaultSC, defaulting to {} if none defined) and standard Hamiltonian defaults (alpha=1, V=0, g=Function[z,-1/z]). Additional optional arguments override each default in order: sc, alpha, V, g. Pass sc={} explicitly to force no switching costs. entries={{vertex,flow},...}, exits={{vertex,cost},...}, sc={{i,k,j,cost},...}. Returns $Failed for unknown keys.

## graphScenario

graphScenario[graph, entries, exits] creates a scenario from any WL directed Graph object. Vertex labels must be positive integers; non-integer labels are not accepted (makeScenario returns Failure). Optional: sc, alpha, V, g.

## gridScenario

gridScenario[dims, entries, exits] creates a scenario on a directed GridGraph[dims]. {n} gives a chain with vertices 1..n; {r,c} gives an r×c grid with vertices 1..r*c (row-major). Optional: sc (switching costs, default {}), alpha, V, g (Hamiltonian defaults from $DefaultHamiltonian).

## makeUnknowns

makeUnknowns[s] returns unknowns[<|"Js" -> ..., "Jts" -> ..., "Us" -> ...|>] for scenario s. "Js" are flow unknowns j[v,e], "Jts" are transition-flow unknowns j[v,e1,e2], and "Us" are value-function unknowns u[v,e].

## unknowns

unknowns[assoc] is the typed head for MFGraphs unknown-variable bundles.

## unknownsData

unknownsData[u, key] returns key from unknowns object u; unknownsData[u] returns the underlying association.

## unknownsQ

unknownsQ[x] returns True iff x is an unknowns[assoc_Association] object.

## altFlowOp

altFlowOp[j][list] returns the alternative: j@@list ==0 || j@@Reverse@list ==0.

## altSwitch

altSwitch[j, u, switchingCosts][v, e1, e2] returns the complementarity condition:
(j[v, e1, e2] == 0) || (u[v, e1] == u[e2, e1] + switchingCosts[{v, e1, e2}])

## consistentSwitchingCosts

consistentSwitchingCosts[switchingcosts][{a,b,c}->S]
returns True if S, the cost of switching from the edge ab to cb, is smaller than any other combination,
such as, ab to bd and then from db to bc.
Returns the condition for this switching cost to satisfy the triangle inequality when S, and the other
switching costs too, does not have a numerical value.

## flowGathering

flowGathering[auxTriples_List][x_] returns the triples that end with x.

## flowSplitting

flowSplitting[AT][UndirectedEdge[a, b]] returns the splitting that start with {a,b}.

## getKirchhoffLinearSystem

getKirchhoffLinearSystem[sys] returns the entry current vector, Kirchhoff matrix, and the variables in the order corresponding to the Kirchhoff matrix.

## getKirchhoffMatrix

getKirchhoffMatrix[sys] returns the entry current vector, Kirchhoff matrix, (critical congestion) cost function placeholder, and the variables in the order corresponding to the Kirchhoff matrix. The third slot is retained for backward compatibility.

## ineqSwitch

ineqSwitch[u, switchingCosts][v, e1, e2] returns the optimality condition at the vertex v related to switching from e1 to e2. Namely,
u[v, e1] <= u[e2, e1] + switchingCosts[{v, e1, e2}]

## isSwitchingCostConsistent

isSwitchingCostConsistent[List of switching costs] is True if all switching costs satisfy the triangle inequality. If some switching costs are symbolic, then it returns the consistency conditions.

## makeSystem

makeSystem[s_scenario, unk_unknowns] constructs an mfgSystem by building the structural equations (SignedFlows, Balance equations, HJ conditions, etc.) from the provided scenario and symbolic unknowns. makeSystem[s_scenario] automatically derives unknowns using makeUnknowns[s].

## mfgBoundaryData

mfgBoundaryData[assoc] is a typed record for boundary conditions and rules.

## mfgBoundaryDataQ

mfgBoundaryDataQ[x] returns True if x is a typed mfgBoundaryData object.

## mfgComplementarityData

mfgComplementarityData[assoc] is a typed record for complementarity conditions and switching costs.

## mfgComplementarityDataQ

mfgComplementarityDataQ[x] returns True if x is a typed mfgComplementarityData object.

## mfgFlowData

mfgFlowData[assoc] is a typed record for flow balance and non-negativity constraints.

## mfgFlowDataQ

mfgFlowDataQ[x] returns True if x is a typed mfgFlowData object.

## mfgHamiltonianData

mfgHamiltonianData[assoc] is a typed record for Hamiltonian residuals and general equations.

## mfgHamiltonianDataQ

mfgHamiltonianDataQ[x] returns True if x is a typed mfgHamiltonianData object.

## mfgSystem

mfgSystem[assoc] is the typed head for an MFG structural equation system. Use makeSystem to construct one; use systemData to access keys.

## mfgSystemQ

mfgSystemQ[x] returns True if x is a typed mfgSystem[assoc_Association] object, False otherwise.

## roundValues

roundValues[x] rounds numerical values in x to a standard precision (10^-10).

## systemData

systemData[sys, key] returns the value associated with key in the system sys, or Missing["KeyAbsent", key] if absent. systemData[sys] returns the underlying Association.

## systemDataFlatten

systemDataFlatten[sys] returns a single flat Association containing all keys from all nested typed sub-records within the system. Useful for backward compatibility with legacy solvers.

## booleanReduceSystem

booleanReduceSystem[sys] solves the mfgSystem sys by converting the preprocessed constraint system to DNF via BooleanConvert, then calling Reduce independently on each disjunct. Each disjunct is a pure conjunction (no Or), so Reduce avoids case-splitting. Non-False results are collected; if the system has a unique equilibrium all non-False results are equivalent. Returns a list of rules when fully determined, or <|"Rules" -> rules, "Equations" -> residual|> when underdetermined. Fails for non-critical congestion systems where Alpha != 1 on any edge. Options: "DisjunctTimeout" (default 30s per Reduce call), "ReturnAll" (default False; True returns all non-False parsed results).

## dnfReduce

dnfReduce[xp, sys] simplifies xp && sys by solving equalities, substituting their solutions throughout the system, and distributing over disjunctions. Returns a DNF expression with all equalities eliminated where possible. dnfReduce[xp, sys, elem] is the 3-argument form used internally to process one conjunct elem from sys.

## dnfReduceSystem

dnfReduceSystem[sys] solves the mfgSystem sys using linear preprocessing followed by dnfReduce instead of Reduce. Handles cases where Reduce times out by using equality-substitution and disjunction-distribution. Returns a list of rules when fully determined, or <|"Rules" -> rules, "Equations" -> residual|> when underdetermined. Fails for non-critical congestion systems where Alpha != 1 on any edge.

## isValidSystemSolution

isValidSystemSolution[sys, sol] checks whether sol (the output of reduceSystem[sys]) satisfies the constraint blocks of sys. Returns True or False. With option "ReturnReport" -> True, returns a detailed association with per-block results. Tolerance for numeric checks is set via "Tolerance" (default 10^-6). For underdetermined solutions the partial rules are checked; blocks that remain symbolic after substitution are reported as Indeterminate, not False.

## reduceSystem

reduceSystem[sys] reduces the structural equations, flow balance, non-negativity constraints, and complementarity conditions of the mfgSystem sys using Reduce over the Reals. Includes AltOptCond (switching-cost optimality complementarity) and IneqSwitchingByVertex (switching-cost optimality inequalities). Fails for non-critical congestion systems where Alpha != 1 on any edge. Returns a list of rules when the system is fully determined, or <|"Rules" -> rules, "Equations" -> residual|> when underdetermined.

## augmentAuxiliaryGraph

augmentAuxiliaryGraph[sys] constructs the road-traffic augmented infrastructure graph from a system's AuxPairs and AuxTriples. Flow edges correspond to j[a,b] and run from {b,a} to {a,b}; transition edges correspond to j[r,i,w] and run from {r,i} to {w,i}. Returns an Association with Graph, Vertices, FlowEdges, TransitionEdges, EdgeVariables, and EdgeKinds.

## mfgAugmentedPlot

mfgAugmentedPlot[s, sys, sol, opts] plots the augmented infrastructure graph (Paper scheme) built by augmentAuxiliaryGraph. Nodes represent road-traffic edge-vertex states, flow edges correspond to j[a,b], and transition edges correspond to j[r,i,w]. Value functions (u) are vertex labels where available, and flow/transition values (j) are edge labels. Supports PlotLabel, GraphLayout, and ImageSize options; the legacy fourth positional title is still accepted.

## mfgFlowPlot

mfgFlowPlot[s, sys, sol, opts] plots a flow-only solution graph with real and auxiliary edges. Edges are displayed as directed, and edge labels show only j flow values. Supports PlotLabel, GraphLayout, and ImageSize options; the legacy fourth positional title is still accepted.

## mfgSolutionPlot

mfgSolutionPlot[s, sys, sol, opts] plots a combined solution graph with real and auxiliary edges. Edge labels include both j and u values, and real-edge directions follow solved net flow orientation. Supports PlotLabel, GraphLayout, and ImageSize options; the legacy fourth positional title is still accepted.

## mfgTransitionPlot

mfgTransitionPlot[s, sys, sol, opts] plots the transition graph of the solution. Nodes are AuxPair states, and directed edges represent transition flows j[r,i,w] from {r,i} to {i,w}. Nodes are labeled with internal values u where available. Supports PlotLabel, GraphLayout, and ImageSize options; the legacy fourth positional title is still accepted.

## scenarioTopologyPlot

scenarioTopologyPlot[s, sys, opts] plots the scenario topology using vertex coloring for entry, exit, and internal vertices. Supports PlotLabel, GraphLayout, and ImageSize options; the legacy third positional title is still accepted.
