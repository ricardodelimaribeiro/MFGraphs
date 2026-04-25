# MFGraphs API Reference

This document is automatically generated from the Wolfram Language package usage metadata.

## alpha

alpha[edge] is the congestion exponent for an edge. Default is 1.

## AltFlowOp

AltFlowOp[j][list] returns the alternative: j@@list ==0 || j@@Reverse@list ==0.

## AltSwitch

AltSwitch[j, u, switchingCosts][v, e1, e2] returns the complementarity condition:
(j[v, e1, e2] == 0) || (u[e1, v] == u[e2, v] + switchingCosts[{v, e1, e2}])
Note that u[r, i] represents the value-function at vertex i on the edge e_{r,i}; thus u[e1, v] and u[e2, v] are both evaluated at vertex v.

## AMScenario

AMScenario[vl, am, entries, exits] creates a scenario from an explicit vertices list vl and adjacency matrix am. Optional: sc, alpha, V, g.

## BuildAuxiliaryTopology

BuildAuxiliaryTopology[model] returns an association with the auxiliary graph and metadata derived from the raw model.

## BuildAuxTriples

BuildAuxTriples[auxGraph] returns the list of all possible {v_in, v_mid, v_out} transitions (triples) in the graph.

## completeScenario

completeScenario[s] fills in derived fields and supplies default Benchmark values ("Tier" -> "core", "Timeout" -> 300) when missing. Returns a new scenario object.

## ConsistentSwitchingCosts

ConsistentSwitchingCosts[switchingcosts][{a,b,c}->S]
returns True if S, the cost of switching from the edge ab to cb, is smaller than any other combination,
such as, ab to bd and then from db to bc.
Returns the condition for this switching cost to satisfy the triangle inequality when S, and the other
switching costs too, does not have a numerical value.

## Cost

Cost[m, edge] is the congestion cost function.

## CycleScenario

CycleScenario[n, entries, exits] creates a scenario on a directed n-cycle (1->2->...->n->1), vertices 1..n. Optional: sc, alpha, V, g.

## DeriveAuxPairs

DeriveAuxPairs[topology] returns the list of all directed edge pairs {u, v} in the auxiliary graph (including entry/exit and reversed graph edges).

## EnsureParallelKernels

EnsureParallelKernels[] launches parallel subkernels if none are running.

## FlowGathering

FlowGathering[auxTriples_List][x_] returns the triples that end with x.

## FlowSplitting

FlowSplitting[AT][UndirectedEdge[a, b]] returns the splitting that start with {a,b}.

## GetExampleScenario

GetExampleScenario[n] returns a 6-arg factory Function[{entries,exits,sc,alpha,V,g}, scenario[...]] for built-in example n. Topology is baked in; all parameters are caller-supplied. GetExampleScenario[n, entries, exits] calls the factory using the canonical switching costs for that case (sc=Automatic resolves via $CaseDefaultSC, defaulting to {} if none defined) and standard Hamiltonian defaults (alpha=1, V=0, g=Function[z,-1/z]). Additional optional arguments override each default in order: sc, alpha, V, g. Pass sc={} explicitly to force no switching costs. entries={{vertex,flow},...}, exits={{vertex,cost},...}, sc={{i,k,j,cost},...}. Returns $Failed for unknown keys.

## GetKirchhoffLinearSystem

GetKirchhoffLinearSystem[sys] returns the entry current vector, Kirchhoff matrix, and the variables in the order corresponding to the Kirchhoff matrix.

## GetKirchhoffMatrix

GetKirchhoffMatrix[sys] returns the entry current vector, Kirchhoff matrix, (critical congestion) cost function placeholder, and the variables in the order corresponding to the Kirchhoff matrix. The third slot is retained for backward compatibility.

## GraphScenario

GraphScenario[graph, entries, exits] creates a scenario from any WL directed Graph object. Optional: sc, alpha, V, g.

## GridScenario

GridScenario[dims, entries, exits] creates a scenario on a directed GridGraph[dims]. {n} gives a chain with vertices 1..n; {r,c} gives an r×c grid with vertices 1..r*c (row-major). Optional: sc (switching costs, default {}), alpha, V, g (Hamiltonian defaults from $DefaultHamiltonian).

## IneqSwitch

IneqSwitch[u, switchingCosts][v, e1, e2] returns the optimality condition at the vertex v related to switching from e1 to e2. Namely,
u[e1, v] <= u[e2, v] + switchingCosts[{v, e1, e2}]
Note that u[r, i] represents the value-function at vertex i on the edge e_{r,i}; thus u[e1, v] and u[e2, v] are both evaluated at vertex v.

## IsSwitchingCostConsistent

IsSwitchingCostConsistent[List of switching costs] is True if all switching costs satisfy the triangle inequality. If some switching costs are symbolic, then it returns the consistency conditions.

## j

j[v, e] or j[v, e1, e2] represents a flow variable.

## makeScenario

makeScenario[assoc] constructs a typed scenario from a raw association. The input must contain a "Model" key whose value is a network topology association accepted by the core scenario/system kernels (required keys: "Vertices List", "Adjacency Matrix", "Entrance Vertices and Flows", "Exit Vertices and Terminal Costs", "Switching Costs"). If "Model" contains a Wolfram Graph under key "Graph", missing "Vertices List" and/or "Adjacency Matrix" are derived automatically. Optional keys: "Hamiltonian" (<|"Alpha" -> a, "V" -> v, "G" -> g, "EdgeAlpha" -> <|{u,v} -> a_uv, ...|>, "EdgeV" -> <|{u,v} -> v_uv, ...|>, "EdgeG" -> <|{u,v} -> g_uv, ...|>|>), "Identity" (name, version), "Benchmark" (tier, timeout), "Visualization", "Inheritance". Default Hamiltonian is Alpha=1 and V=0 on all edges, with G[z]=-1/z (overridable globally and per edge). Boundary and switching-cost values must be numeric. Returns a scenario[...] object on success or Failure[...] on error.

## makeSystem

makeSystem[s_scenario, unk_unknowns] constructs an mfgSystem by building the structural equations (SignedFlows, Balance equations, HJ conditions, etc.) from the provided scenario and symbolic unknowns. makeSystem[s_scenario] automatically derives unknowns using makeUnknowns[s].

## makeUnknowns

makeUnknowns[s] returns unknowns[<|"js" -> ..., "jts" -> ..., "us" -> ...|>] for scenario s. "js" are flow unknowns j[v,e], "jts" are transition-flow unknowns j[v,e1,e2], and "us" are value-function unknowns u[v,e].

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

## MFGParallelMap

MFGParallelMap[f, list] applies f to each element of list, using ParallelMap
when Length[list] >= $MFGraphsParallelThreshold (and kernels are already running)
or Length[list] >= $MFGraphsParallelLaunchThreshold (when kernels must be launched).
This avoids paying ~3s kernel launch overhead for small workloads.

## MFGPrint

MFGPrint[args___] prints args only when $MFGraphsVerbose is True.

## MFGPrintTemporary

MFGPrintTemporary[args___] prints a temporary message only when $MFGraphsVerbose is True.

## mfgSystem

mfgSystem[assoc] is the typed head for an MFG structural equation system. Use makeSystem to construct one; use SystemData to access keys.

## mfgSystemQ

mfgSystemQ[x] returns True if x is a typed mfgSystem[assoc_Association] object, False otherwise.

## RoundValues

RoundValues[x] rounds numerical values in x to a standard precision (10^-10).

## scenario

scenario[assoc] is the typed head for a MFGraphs scenario object. Use makeScenario to construct one; use ScenarioData to access keys.

## ScenarioData

ScenarioData[s, key] returns the value associated with key in the scenario s, or Missing["KeyAbsent", key] if absent. ScenarioData[s] returns the underlying Association.

## scenarioQ

scenarioQ[x] returns True if x is a typed scenario[assoc_Association] object, False otherwise.

## SystemData

SystemData[sys, key] returns the value associated with key in the system sys, or Missing["KeyAbsent", key] if absent. SystemData[sys] returns the underlying Association.

## SystemDataFlatten

SystemDataFlatten[sys] returns a single flat Association containing all keys from all nested typed sub-records within the system. Useful for backward compatibility with legacy solvers.

## u

u[v, e] represents a utility/potential variable. u[r, i] represents the value of the value-function at vertex v_i on the edge e_{r,i}.

## unknowns

unknowns[assoc] is the typed head for MFGraphs unknown-variable bundles.

## UnknownsData

UnknownsData[u, key] returns key from unknowns object u; UnknownsData[u] returns the underlying association.

## unknownsQ

unknownsQ[x] returns True iff x is an unknowns[assoc_Association] object.

## validateScenario

validateScenario[s] checks that the scenario s has all required Model keys and that the Model value is an Association. Returns s unchanged on success, or Failure["ScenarioValidation", <|"Message" -> msg, "MissingKeys" -> {...}|>] on failure.

## z

z[v] represents a vertex potential variable.

## $MFGraphsParallelReady

False is the symbol for the Boolean value false. 

## $MFGraphsVerbose

False is the symbol for the Boolean value false. 
