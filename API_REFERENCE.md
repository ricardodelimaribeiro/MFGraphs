# MFGraphs API Reference

This document is automatically generated from the Wolfram Language package usage metadata.

## alpha

alpha[edge] is the congestion exponent for an edge. Default is 1.

## Cost

Cost[m, edge] is the congestion cost function.

## ensureParallelKernels

ensureParallelKernels[] launches parallel subkernels if none are running.

## j

j[a, b] is the flow on edge {a,b} from a to b. j[r, i, w] is the fraction of the flow j[r,i] that transitions to edge e_{i,w} at junction i.

## mfgPrint

mfgPrint[args___] prints args only when $MFGraphsVerbose is True.

## mfgPrintTemporary

mfgPrintTemporary[args___] prints a temporary message only when $MFGraphsVerbose is True.

## u

u[a, b] is the value of the value function at vertex b of edge {a,b}.

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

makeScenario[assoc] constructs a typed scenario from a raw association. The input must contain a "Model" key whose value is a network topology association (required keys: "Vertices", "Adjacency", "Entries", "Exits", "Switching"). If "Model" contains a Wolfram Graph under key "Graph", missing "Vertices" and/or "Adjacency" are derived automatically. Optional keys: "Hamiltonian" (<|"Alpha" -> a, "V" -> v, "G" -> g, "EdgeAlpha" -> <|{u,v} -> a_uv, ...|>, "EdgeV" -> <|{u,v} -> v_uv, ...|>, "EdgeG" -> <|{u,v} -> g_uv, ...|>|>), "Identity" (name, version), "Benchmark" (tier, timeout), "Visualization", "Inheritance". Default Hamiltonian is Alpha=1 and V=-1 on all edges, with G[z]=-1/z (overridable globally and per edge). V/G/EdgeV/EdgeG are validated and preserved for future density and visualization work, but current system construction applies only Alpha/EdgeAlpha. Boundary values must be numeric; switching-cost values must be numeric or Infinity. Returns a scenario[...] object on success or Failure[...] on error.

## scenario

scenario[assoc] is the typed head for a MFGraphs scenario object. Use makeScenario to construct one; use scenarioData to access keys.

## scenarioData

scenarioData[s, key] returns the value associated with key in the scenario s, or Missing["KeyAbsent", key] if absent. scenarioData[s] returns the underlying Association.

## scenarioQ

scenarioQ[x] returns True if x is a typed scenario[assoc_Association] object, False otherwise.

## validateScenario

validateScenario[s] checks that the scenario s has all required Model keys and that the Model value is an Association. Returns s unchanged on success, or Failure["ScenarioValidation", <|"Message" -> msg, "MissingKeys" -> {...}|>] on failure.

## amScenario

amScenario[vl, am, entries, exits] creates a scenario from an explicit vertices list vl and adjacency matrix am. Vertex labels in vl must be positive integers; non-integer labels are not accepted. The adjacency matrix is treated as network connections and symmetrized before system construction. Optional: sc, alpha, V, g.

## cycleScenario

cycleScenario[n, entries, exits] creates a scenario on n-cycle connections, vertices 1..n. Topology is symmetrized into undirected network edges before system construction. Optional: sc, alpha, V, g.

## getExampleScenario

getExampleScenario[n] returns a 6-arg factory Function[{entries,exits,sc,alpha,V,g}, scenario[...]] for built-in example n (integer 1-23 or named string). Topology is baked in; all parameters are caller-supplied. getExampleScenario[n, entries, exits] calls the factory with canonical defaults:   sc=Automatic resolves via $CaseDefaultSC (falls back to {}),   alpha=1 (critical congestion), V=0, g=0 (Hamiltonian parameters passed through but not yet applied in system construction). getExampleScenario[n, entries, exits, sc] overrides switching costs (pass {} for none). getExampleScenario[n, entries, exits, sc, alpha] also overrides alpha. getExampleScenario[n, entries, exits, sc, alpha, V] also overrides V. getExampleScenario[n, entries, exits, sc, alpha, V, g] overrides all parameters. entries={{vertex,flow},...}; exits={{vertex,cost},...}; sc={{i,k,j,cost},...} or Association with 3-tuple keys. Returns $Failed for unknown keys.

## getExampleScenarioMetadata

getExampleScenarioMetadata[key] returns metadata for a built-in example scenario, including source/provenance details when available. Returns $Failed for unknown keys.

## graphScenario

graphScenario[graph, entries, exits] creates a scenario from any WL Graph object. Graph edges are treated as network connections and symmetrized before system construction; movement restrictions should be encoded with switching costs such as Infinity. Vertex labels must be positive integers; non-integer labels are not accepted (makeScenario returns Failure). Optional: sc, alpha, V, g.

## gridScenario

gridScenario[dims, entries, exits] creates a scenario on GridGraph[dims] connections. {n} gives a chain with vertices 1..n; {r,c} gives an r×c grid with vertices 1..r*c (row-major). Topology is symmetrized into undirected network edges before system construction. Optional: sc (switching costs, default {}), alpha, V, g (Hamiltonian defaults are supplied by makeScenario; V/G are preserved but not applied by current system construction).

## makeSymbolicUnknowns

makeSymbolicUnknowns[s] returns symbolicUnknowns[<|"Js" -> ..., "Jts" -> ..., "Us" -> ...|>] for scenario s. "Js" are exact flow variables j[a,b], "Jts" are exact transition-flow variables j[r,i,w], and "Us" are exact value-function variables u[a,b].

## symbolicUnknowns

symbolicUnknowns[assoc] is the typed head for topology-derived exact symbolic variable bundles used by MFGraphs structural graph systems.

## symbolicUnknownsData

symbolicUnknownsData[u, key] returns key from symbolicUnknowns object u; symbolicUnknownsData[u] returns the underlying association.

## symbolicUnknownsQ

symbolicUnknownsQ[x] returns True iff x is a symbolicUnknowns[assoc_Association] object.

## unknown

unknown is reserved for future numeric MFGraphs solver fields over grids or layouts; it is not used by the exact symbolic graph-system constructors.

## unknowns

unknowns is reserved for future numeric collections of unknown fields over grids or layouts, following the Maydan-style numeric solver boundary. Use symbolicUnknowns for current exact symbolic graph systems.

## altFlowOp

altFlowOp[j][list] returns the alternative: j@@list ==0 || j@@Reverse@list ==0.

## altSwitch

altSwitch[j, u, switchingCosts][r, i, w] returns the complementarity condition at junction i for the transition from r to w:
(j[r, i, w] == 0) || (u[w, i] + switchingCosts[r, i, w] - u[r, i] == 0).

## buildBoundaryData

buildBoundaryData[s, topology] builds typed boundary equations, entry rules, exit inequalities, and entry/exit metadata. IneqExitValues stores upper-bound inequalities u[auxExit,N]<=cost for each exit; AltExitCond stores the complementarity conditions j[N,auxExit]==0||u[auxExit,N]==cost.

## buildComplementarityData

buildComplementarityData[s, topology, unk] builds typed complementarity alternatives and switching inequalities.

## buildFlowData

buildFlowData[s, topology, unk] builds typed flow-balance equations and non-negativity constraints.

## buildHamiltonianData

buildHamiltonianData[s, topology, flowData] builds typed Hamiltonian residual equations for the system. EqGeneral encodes the edge-level HJB equation u[a,b]-u[b,a]+j[a,b]-j[b,a] = nrhs unconditionally for every undirected edge {a,b}, where nrhs = 0 for Alpha==1 and m - Sign[m] m^alpha otherwise. The equation is enforced even when net flow is zero (zero-flow edges force u[a,b]=u[b,a]), which propagates through switching inequalities to pin value variables at bypassed exit nodes to their terminal cost. Current system construction uses Alpha/EdgeAlpha; V/G/EdgeV/EdgeG are preserved on scenarios for future work.

## flowGathering

flowGathering[auxTriples_List][x_] returns the triples that end with x.

## flowSplitting

flowSplitting[AT][UndirectedEdge[a, b]] returns the splitting that start with {a,b}.

## getKirchhoffLinearSystem

getKirchhoffLinearSystem[sys] returns the entry current vector, Kirchhoff matrix, and the variables in the order corresponding to the Kirchhoff matrix.

## getKirchhoffMatrix

getKirchhoffMatrix[sys] returns the entry current vector, Kirchhoff matrix, (critical congestion) cost function placeholder, and the variables in the order corresponding to the Kirchhoff matrix. The third slot is retained for backward compatibility.

## ineqSwitch

ineqSwitch[u, switchingCosts][r, i, w] returns the optimality condition at junction i for the transition from r to w. Namely,
u[w, i] + switchingCosts[r, i, w] - u[r, i] >= 0.

## makeSystem

makeSystem[s_scenario, unk_symbolicUnknowns] constructs an mfgSystem by building the structural equations (SignedFlows, Balance equations, HJ conditions, etc.) from the provided scenario and exact symbolic unknown bundle. makeSystem[s_scenario] automatically derives symbolicUnknowns using makeSymbolicUnknowns[s].

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

## switchingCostLookup

switchingCostLookup[sc] returns a function f such that f[r, i, w] gives the switching cost for the transition from e_{r,i} to e_{i,w}, defaulting to 0 if absent.

## systemData

systemData[sys, key] returns the value associated with key in the system sys, or Missing["KeyAbsent", key] if absent. systemData[sys] returns the underlying Association.

## systemDataFlatten

systemDataFlatten[sys] returns a single flat Association containing all keys from all nested typed sub-records within the system. Useful for backward compatibility with legacy solvers.

## activeSetReduceSystem

activeSetReduceSystem[sys] is an opt-in exact active-set solver for the critical-congestion linear complementarity structure. It enumerates small complementarity alternatives incrementally with exact linear substitution and falls back to the proven exact DNF reducer for larger residual variable sets. Returns the same rule/residual shape as dnfReduceSystem. Fails for non-critical congestion systems where Alpha != 1 on any edge.

## booleanMinimizeReduceSystem

booleanMinimizeReduceSystem[sys] solves the mfgSystem sys by attacking the disjunctive structure of the preprocessed constraint system before DNF expansion. It (1) prunes individual complementarity arms that are infeasible against the linear part via FindInstance; (2) decomposes the surviving disjunctive atoms into connected components by shared variables; (3) BooleanMinimizes each component to minimal DNF and Reduces per disjunct. Returns the same rule/residual shape as booleanReduceSystem. Fails for non-critical congestion systems where Alpha != 1 on any edge. Options: "ArmTimeout" (default 2s per FindInstance arm check), "DisjunctTimeout" (default 30s per Reduce call), "ReturnAll" (default False).

## booleanMinimizeSystem

booleanMinimizeSystem[sys] is a head-to-head variant of booleanReduceSystem that calls BooleanMinimize[constraints, "DNF"] in place of BooleanConvert[constraints, "DNF"]. This is exact minimal-DNF Boolean minimization, analogous to classical two-level SOP minimization such as Quine-McCluskey/Petrick-style methods, followed by Reduce per disjunct. Wolfram does not document BooleanMinimize as a specific QMC/Petrick implementation. Same preprocessing and return shape as booleanReduceSystem. Use to compare the two Boolean-stage operations on identical input. Fails for non-critical congestion systems where Alpha != 1 on any edge. Options: "DisjunctTimeout" (default 30s per Reduce call), "ReturnAll" (default False).

## booleanReduceSystem

booleanReduceSystem[sys] solves the mfgSystem sys by converting the preprocessed constraint system to DNF via BooleanConvert, then calling Reduce independently on each disjunct. This is DNF conversion followed by real quantifier elimination / CAD-style solving per pure conjunction. Each disjunct has no Or, so Reduce avoids case-splitting. Non-False results are collected; if the system has a unique equilibrium all non-False results are equivalent. Returns a list of rules when fully determined, or <|"Rules" -> rules, "Residual" -> residual|> when underdetermined. Fails for non-critical congestion systems where Alpha != 1 on any edge. Options: "DisjunctTimeout" (default 30s per Reduce call), "ReturnAll" (default False; True returns all non-False parsed results).

## computeKirchhoffResidual

computeKirchhoffResidual[sys, sol] calculates the Kirchhoff Residual (the maximum absolute divergence in the transition flux conservation equations) for the solution sol. Returns the maximum numerical residual, or Missing["NotAvailable"] if the solution is symbolic or incomplete.

## directCriticalSystem

directCriticalSystem[sys] is an explicit opt-in solver for critical congestion systems with all numeric zero switching costs and equal numeric exit costs. It uses graph distance to exits to forbid flow directions that move farther from every exit, then solves the resulting feasibility system. Returns the standard raw solver shape or Failure["directCriticalSystem", ...].

## dnfReduce

dnfReduce[xp, sys] simplifies xp && sys by solving equalities, substituting their solutions throughout the system, and distributing over disjunctions. Returns a DNF expression with all equalities eliminated where possible. dnfReduce[xp, sys, elem] is the 3-argument form used internally to process one conjunct elem from sys. Implemented with direct recursion: the Or case spawns one recursive call per branch and the Equal case recurses on the substituted remainder. Deeply nested Or-chains may approach $RecursionLimit; see dnfReduceProcedural for a stack-based iterative alternative.

## dnfReduceSystem

dnfReduceSystem[sys] solves the mfgSystem sys using linear preprocessing followed by dnfReduce instead of Reduce. Handles cases where Reduce times out by using equality-substitution and disjunction-distribution. Returns a list of rules when fully determined, or <|"Rules" -> rules, "Residual" -> residual|> when underdetermined. Fails for non-critical congestion systems where Alpha != 1 on any edge.

## findInstanceSystem

findInstanceSystem[sys] solves the mfgSystem sys by collecting and linearly preprocessing constraints, then calling FindInstance over the remaining real variables. This is real satisfiability / instance finding using Wolfram's real-system solver backend. Returns one feasible list of rules. If no instance is found or the final solve times out, returns <|"Rules" -> accumulatedRules, "Residual" -> False|>. Fails for non-critical congestion systems where Alpha != 1 on any edge. Options: "Timeout" (default Infinity).

## flowFirstCriticalSystem

flowFirstCriticalSystem[sys] is an explicit opt-in solver for critical congestion systems. It minimizes the critical quadratic flow objective over Kirchhoff/nonnegative flow constraints, then recovers value variables from the remaining system. Returns the standard raw solver shape or Failure["flowFirstCriticalSystem", ...].

## isValidSystemSolution

isValidSystemSolution[sys, sol] checks whether sol (the output of reduceSystem[sys]) satisfies the constraint blocks of sys. Returns True or False. With option "ReturnReport" -> True, returns a detailed association with per-block results. Tolerance for numeric checks is set via "Tolerance" (default 10^-6). For underdetermined solutions the partial rules are checked; blocks that remain symbolic after substitution are reported as Indeterminate, not False.

## optimizedDNFReduceSystem

optimizedDNFReduceSystem[sys] is an opt-in exact solver for critical congestion systems. It follows the same preprocessing and output contract as dnfReduceSystem. It carries small DNF branch families directly as rules plus residual constraints, and falls back to the proven exact DNF reducer for larger residual variable sets. Returns a list of rules when fully determined, or <|"Rules" -> rules, "Residual" -> residual|> when underdetermined. Fails for non-critical congestion systems where Alpha != 1 on any edge.

## reduceSystem

reduceSystem[sys] reduces the structural equations, flow balance, non-negativity constraints, and complementarity conditions of the mfgSystem sys using Wolfram Reduce over the Reals; for polynomial real systems this is the CAD/real quantifier-elimination backend. Includes AltOptCond (switching-cost optimality complementarity) and IneqSwitchingByVertex (switching-cost optimality inequalities). Fails for non-critical congestion systems where Alpha != 1 on any edge. Returns a list of rules when the system is fully determined, or <|"Rules" -> rules, "Residual" -> residual|> when underdetermined.

## solutionBranchCostReport

solutionBranchCostReport[sys, sol] ranks top-level residual branches by a diagnostic critical-congestion objective. It reports terminal exit cost, switching cost, edge-flow quadratic cost, total objective, determined and residual variables, and validation status for each branch.

## solutionReport

solutionReport[sys, sol] returns a read-only diagnostic association for an existing solver result. The report includes result kind, validation report, Kirchhoff residual, residual branch diagnostics, primary residual variables, and transition-flow determinacy diagnostics. It does not rewrite sol.

## clearSolveCache

clearSolveCache[] empties solveScenario's session-scoped memoization cache. Call between benchmark passes to measure cold-start cost, or to release memory in long-running sessions.

## SolveMFG

SolveMFG[s] solves a typed scenario object by delegating to solveScenario. SolveMFG[assoc] provides backward compatibility for legacy raw-association solving. It constructs a scenario and delegates to solveScenario. SolveMFG[s, solver] or SolveMFG[assoc, solver] uses the specified solver function.

## solveScenario

solveScenario[s] automatically constructs exact symbolic unknowns, builds the structural system, and calls dnfReduceSystem. solveScenario[{s1, s2, ...}] solves multiple populations (scenarios) and returns a list of solutions. solveScenario[..., solver] uses the specified solver function (e.g., reduceSystem). Results are memoized per (scenario, solver) within the session; call clearSolveCache[] to wipe.

## augmentAuxiliaryGraph

augmentAuxiliaryGraph[sys] constructs the road-traffic augmented infrastructure graph from a system's AuxPairs and AuxTriples. Returns an Association containing the Graph, Vertices, FlowEdges, TransitionEdges, EdgeVariables, and EdgeKinds.

## BendFactor

BendFactor is an option for richNetworkPlot.

## rawNetworkPlot

rawNetworkPlot[s, sys, opts] and rawNetworkPlot[s, sys, sol, opts] render the physical network. Real network vertices are gray. Auxiliary entry/exit vertices and boundary edges are hidden by default. Options: ShowAuxiliaryVertices (default False), ShowBoundaryData (default False; when True forces ShowAuxiliaryVertices), ShowBoundaryValues (default True; shows u-values on auxiliary exit vertices when boundary data is shown), ShowFlowLabels (default Automatic; True when sol provided), ShowValueLabels (default Automatic; True when sol provided; shows u-values with 1/4, 1/2, 3/4 interpolations), ShowDensityLabels (default False; shows inferred density m), ColorFunction (default Automatic; RedBlue blend over u-values), ShowLegend (default True), GraphLayout (default Automatic), PlotLabel (default Automatic), ImageSize (default Large). With a solution provided, auxiliary vertices (when shown) and shown vertex states are colored by u-value gradient; physical entry/exit vertices remain gray.

## richNetworkPlot

richNetworkPlot[s, sys, opts] and richNetworkPlot[s, sys, sol, opts] render the augmented state-space graph. Nodes are pairs {a,b} representing oriented edge states; edges are flow arcs j[a,b] (blue) and transition arcs j[r,i,w] (red), drawn as quadratic Bezier curves so anti-parallel pairs separate. Only nodes whose label position (b) is an auxiliary entry/exit vertex receive boundary colors. Options: ShowFlowEdges (default True; False yields the transition-only graph), ShowBoundaryData (default False; overlays entry-flow and exit-cost labels), ShowBoundaryValues (default True; shows u/cost on boundary nodes), ShowFlowLabels (default Automatic; True when sol provided), ShowValueLabels (default Automatic; True when sol provided), UseColorFunction (default False; when True colors nodes and eligible edges by u-values using ColorFunction), ColorFunction (default Automatic), ShowLegend (default True), BendFactor (default 0.15; arc curvature as fraction of edge length), GraphLayout (default Automatic), PlotLabel (default Automatic), ImageSize (default Large).

## ShowAuxiliaryVertices

ShowAuxiliaryVertices is an option for rawNetworkPlot.

## ShowBoundaryData

ShowBoundaryData is an option for rawNetworkPlot and richNetworkPlot.

## ShowBoundaryValues

ShowBoundaryValues is an option for rawNetworkPlot and richNetworkPlot.

## ShowDensityLabels

ShowDensityLabels is an option for rawNetworkPlot.

## ShowFlowEdges

ShowFlowEdges is an option for richNetworkPlot.

## ShowFlowLabels

ShowFlowLabels is an option for rawNetworkPlot and richNetworkPlot.

## ShowLegend

ShowLegend is an option for rawNetworkPlot and richNetworkPlot.

## ShowValueLabels

ShowValueLabels is an option for rawNetworkPlot and richNetworkPlot.

## UseColorFunction

UseColorFunction is an option for richNetworkPlot.
