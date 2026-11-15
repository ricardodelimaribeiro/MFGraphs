(* Wolfram Language package *)
(*
   scenarioTools.wl — typed scenario kernel for MFGraphs.

   Provides a typed wrapper around network-plus-parameters data so that
   scenarios can be constructed, validated, completed, and stored uniformly.

   Lifecycle:
     makeScenario[rawAssoc]  →  normalizes, validates, completes, and wraps in scenario[...]
     validateScenario[s]     →  True if s is valid; otherwise issues a message and returns False
     completeScenario[s]     →  fills defaults, caches topology, completes switching costs
     scenarioQ[x]            →  True iff x is a well-formed typed scenario

   Top-level blocks inside scenario[<|...|>]:
     "Model"           — network and boundary data (keys: "Vertices", "Adjacency",
                         "Entries", "Exits", "Switching"; optional "Graph",
                         "EdgeEntryTolls")
     "Hamiltonian"     — "Alpha", "V", "G", "EdgeAlpha", "EdgeV", "EdgeG"
     "Topology"        — derived auxiliary graph data (cached)
     "InputProvenance" — boundary, switching, toll and Hamiltonian input as given
     "Identity"        — name, version
     "Benchmark"       — tier, timeout
     "Visualization", "Inheritance" — optional, passed through unchanged

   Auxiliary vertices: each entry vertex v gets the vertex auxEntry[v] and each
   exit vertex v gets auxExit[v]; boundary flows run along {auxEntry[v], v}
   and {v, auxExit[v]}.

   Switching costs: every admissible transition without an explicit cost gets
   cost 0, then the finite costs at each junction are replaced by their metric
   closure, so c(r,v,w) <= c(r,v,k) + c(k,v,w) holds among them. Blocked
   transitions (Infinity) and flow-dependent costs (pure functions) are left
   unchanged.
*)

BeginPackage["scenarioTools`", {"primitives`", "utilities`"}]

(* Lets a re-evaluation of this file reset the usage messages below. *)
Unprotect[auxEntry, auxExit];

(* --- Public API declarations --- *)

auxEntry::usage =
"auxEntry[v] labels the auxiliary entry vertex attached to entry vertex v. \
Boundary inflow enters along the edge {auxEntry[v], v}. auxEntry is Protected \
and has no definitions; it only marks a vertex's role.";

auxExit::usage =
"auxExit[v] labels the auxiliary exit vertex attached to exit vertex v. \
Boundary outflow leaves along the edge {v, auxExit[v]}. auxExit is Protected \
and has no definitions; it only marks a vertex's role.";

scenario::usage =
"scenario[assoc] is the typed head for a MFGraphs scenario object. Use makeScenario \
to construct one; use scenarioData to access keys.";

scenarioQ::usage =
"scenarioQ[x] returns True if x is a typed scenario[assoc_Association] object, \
False otherwise.";

makeScenario::usage =
"makeScenario[assoc] constructs a typed scenario from a raw association. The input \
must contain a \"Model\" key whose value is a network association with keys \
\"Vertices\" (distinct integers), \"Adjacency\" (square matrix, zero diagonal), \
\"Entries\" ({vertex, inflow} pairs, inflow >= 0), \"Exits\" ({vertex, exit cost} \
pairs), and \"Switching\" ({r, v, w, cost} rows or an Association {r, v, w} -> cost). \
If \"Model\" contains a Wolfram Graph under key \"Graph\", missing \"Vertices\" \
and/or \"Adjacency\" are derived from it. \"Switching\" may be omitted when \
\"EdgeEntryTolls\" (<|{v, w} -> toll, ...|>, a toll for entering edge {v, w} at v) \
is given. Boundary values must be real; each vertex may appear at most once in \
Entries and at most once in Exits, and the two sets must be disjoint. Switching \
costs must be keyed by admissible transitions (consecutive edges, no U-turn) and \
be nonnegative numbers, Infinity (a blocked transition), or pure functions of the \
transition flow. Every admissible transition without an explicit cost or toll gets \
cost 0; the finite costs at each junction are then replaced by their metric closure \
(the cheapest sequence of switches), so they satisfy the triangle inequality \
c(r,v,w) <= c(r,v,k) + c(k,v,w). Blocked transitions (Infinity) and pure-function \
costs are not changed by the closure. \
Optional keys: \"Hamiltonian\" (<|\"Alpha\" -> a, \"V\" -> v, \"G\" -> g, \
\"EdgeAlpha\" -> <|{u,v} -> a_uv, ...|>, \"EdgeV\" -> <|{u,v} -> v_uv, ...|>, \
\"EdgeG\" -> <|{u,v} -> g_uv, ...|>|>; other keys are rejected), \
\"Identity\" (name, version), \"Benchmark\" (tier, timeout), \"Visualization\", \
\"Inheritance\". Default Hamiltonian is Alpha=1 and V=-1 on all edges, with \
G[z]=-1/z (overridable globally and per edge). V/G/EdgeV/EdgeG are validated \
and preserved for density and visualization work; system construction does not \
use them, and the solvers accept only Alpha=1. Returns a scenario[...] object on \
success or Failure[\"ScenarioValidation\", ...] describing the first problem found.";

validateScenario::usage =
"validateScenario[s] returns True if s is a scenario object whose \"Model\" and \
\"Hamiltonian\" pass every check that makeScenario applies, and False otherwise. \
On False it issues validateScenario::invalid (or validateScenario::notscenario) \
naming the first problem found; makeScenario reports the same problems as a Failure.";

completeScenario::usage =
"completeScenario[s] returns s with defaults filled in: absent Hamiltonian values \
take their defaults, Benchmark defaults (\"Tier\" -> \"core\", \"Timeout\" -> 300) \
are added when missing, and switching costs are completed over every admissible \
transition. If s has no cached \"Topology\", it issues completeScenario::notopology, \
rebuilds the topology from the Model, and stores it in the result. completeScenario \
does not validate s; use validateScenario first.";

scenarioData::usage =
"scenarioData[s, key] returns the value associated with key in the scenario s, or \
Missing[\"KeyAbsent\", key] if absent. scenarioData[s] returns the underlying \
Association.";

buildAuxiliaryTopology::usage =
"buildAuxiliaryTopology[model] returns an association with the network graph, the \
auxiliary graph (with auxEntry[v] / auxExit[v] boundary vertices), the directed \
pairs, and the admissible transition triples derived from the model, or $Failed \
if the model's vertices, adjacency, entries, or exits are malformed.";

Begin["`Private`"]

scenarioFailure::usage =
"scenarioFailure[msg] builds a ScenarioValidation Failure. scenarioFailure[msg, details] adds the entries of the details Association.";

undirectedAdjacency::usage =
"undirectedAdjacency[adjacency] returns the 0/1 adjacency matrix of the undirected network underlying adjacency.";

realNumberQ::usage =
"realNumberQ[x] returns True when x is a real numeric quantity; complex values are rejected.";

switchingCostValueQ::usage =
"switchingCostValueQ[c] returns True when c is a nonnegative real number, Infinity (a blocked transition), or a pure Function.";

hamiltonianTermQ::usage =
"hamiltonianTermQ[value] returns True when value is numeric or a pure Function, as allowed for V, G, EdgeV and EdgeG.";

buildAdjacencyFromGraph::usage =
"buildAdjacencyFromGraph[graph, vertices] returns graph's adjacency matrix ordered by vertices, or $Failed for unknown vertices.";

normalizeScenarioModel::usage =
"normalizeScenarioModel[model] fills Vertices and Adjacency from Graph when possible and normalizes the stored Graph; a model that does not match its Graph is returned unchanged for validation to report.";

normalizeSwitchingCosts::usage =
"normalizeSwitchingCosts[sc] converts list switching-cost rows to an Association and leaves Associations unchanged.";

normalizeHamiltonianSpec::usage =
"normalizeHamiltonianSpec[spec] returns spec with the default Hamiltonian values filled in for absent keys.";

modelEdgeSet::usage =
"modelEdgeSet[vertices, adjacency] returns an Association whose keys are the edges {a, b} of the undirected network, in both orientations.";

switchingSpecQ::usage =
"switchingSpecQ[sc] returns True when sc is a list of {r, v, w, cost} rows or an Association {r, v, w} -> cost with valid cost values.";

scenarioInputProblem::usage =
"scenarioInputProblem[assoc] returns a ScenarioValidation Failure for the first problem in assoc's Model or Hamiltonian that can be found without the auxiliary topology, or None.";

hamiltonianProblem::usage =
"hamiltonianProblem[spec, edgeSet] returns a ScenarioValidation Failure for the first problem in a Hamiltonian spec, or None.";

switchingProblem::usage =
"switchingProblem[model, auxTriples] returns a ScenarioValidation Failure listing switching-cost keys that are not admissible transitions, or None.";

checkedTopology::usage =
"checkedTopology[assoc] validates assoc and returns its auxiliary topology (the cached one when present), or the first ScenarioValidation Failure.";

expandEdgeEntryTolls::usage =
"expandEdgeEntryTolls[model] turns each toll on entering edge {v, w} at v into the switching costs {r, v, w} for every other neighbour r of v.";

metricClosureSwitchingCosts::usage =
"metricClosureSwitchingCosts[sc] replaces the finite numeric costs at each junction by their shortest-path closure; Infinity and Function costs are left unchanged.";

junctionClosure::usage =
"junctionClosure[sc, triples] returns rules giving the shortest-path closure of the finite costs sc on the triples that share one junction.";

auxTriplesFromPairs::usage =
"auxTriplesFromPairs[directedPairs] returns the admissible transitions {r, v, w}: consecutive directed pairs {r, v} and {v, w} with r != w.";

completedSwitchingCosts::usage =
"completedSwitchingCosts[model, auxTriples] gives every admissible transition a cost: explicit switching costs, then expanded edge entry tolls, then 0; the finite costs are then replaced by their metric closure.";

(* Required keys that must appear inside the "Model" block. *)
$RequiredModelKeys = {"Vertices", "Adjacency", "Entries", "Exits", "Switching"};

(* Default benchmark settings applied by completeScenario. *)
$DefaultBenchmarkTier    = "core";
$DefaultBenchmarkTimeout = 300;

(* Default Hamiltonian: alpha = 1 and V = -1 on all edges, and g(z) = -1/z.
   Its keys are the only keys a "Hamiltonian" spec may contain. *)
$DefaultHamiltonian = <|
    "Alpha" -> 1,
    "V" -> -1,
    "G" -> Function[z, -1 / z],
    "EdgeAlpha" -> <||>,
    "EdgeV" -> <||>,
    "EdgeG" -> <||>
|>;

(* --- Small predicates and normalizers --- *)

scenarioFailure[msg_String] := scenarioFailure[msg, <||>];
scenarioFailure[msg_String, details_Association] :=
    Failure["ScenarioValidation", Join[<|"Message" -> msg, "MissingKeys" -> {}|>, details]];

undirectedAdjacency[adjacency_] := Unitize[adjacency + Transpose[adjacency]];

realNumberQ[x_] := NumericQ[x] && Im[x] === 0;

(* Infinity is the explicit sentinel for blocked transitions. In WL this is the
   canonical positive infinity value; other directed infinities are intentionally
   not accepted as switching costs. *)
switchingCostValueQ[c_] := (realNumberQ[c] && c >= 0) || c === Infinity || MatchQ[c, _Function];

(* V may be numeric (constant along the edge) or a pure Function of the edge
   coordinate x in [0,1]; G may be numeric or a pure Function of the density.
   System construction uses only Alpha; V and G are consumed downstream for
   density and value-profile recovery (see edgeDensityProfile / edgeValueProfile). *)
hamiltonianTermQ[value_] := NumericQ[value] || MatchQ[value, _Function];

buildAdjacencyFromGraph[graph_Graph, vertices_List] :=
    Module[{index, order},
        index = AssociationThread[VertexList[graph], Range[VertexCount[graph]]];
        order = Lookup[index, vertices, Missing["UnknownVertex"]];
        If[MemberQ[order, _Missing], $Failed, AdjacencyMatrix[graph][[order, order]]]
    ];

normalizeScenarioModel[model_Association] :=
    Module[{graph, vertices, adjacency, n},
        graph = Lookup[model, "Graph"];
        If[!MatchQ[graph, _Graph], Return[model, Module]];
        vertices  = Lookup[model, "Vertices", VertexList[graph]];
        adjacency = Lookup[model, "Adjacency",
            If[ListQ[vertices], buildAdjacencyFromGraph[graph, vertices], $Failed]];
        n = If[ListQ[vertices], Length[vertices], -1];
        (* A model that does not match its Graph is left unchanged, so validation
           reports the problem instead of AdjacencyGraph printing messages. *)
        If[!DuplicateFreeQ[vertices] || !MatrixQ[adjacency] || Dimensions[adjacency] =!= {n, n},
            Return[model, Module]];
        Join[model, <|
            "Vertices"  -> vertices,
            "Adjacency" -> adjacency,
            "Graph"     -> AdjacencyGraph[vertices, undirectedAdjacency[adjacency], DirectedEdges -> False]
        |>]
    ];
normalizeScenarioModel[other_] := other;

normalizeSwitchingCosts[sc_Association] := sc;
normalizeSwitchingCosts[sc_List] := AssociationThread[Most /@ sc, Last /@ sc];

normalizeHamiltonianSpec[spec_Association] := Join[$DefaultHamiltonian, spec];
normalizeHamiltonianSpec[_] := $DefaultHamiltonian;

modelEdgeSet[vertices_List, adjacency_] :=
    AssociationMap[True &,
        Map[vertices[[#]] &, SparseArray[undirectedAdjacency[adjacency]]["NonzeroPositions"]]];

switchingSpecQ[sc_Association] :=
    AllTrue[Keys[sc], MatchQ[{_, _, _}]] && AllTrue[Values[sc], switchingCostValueQ];
switchingSpecQ[sc_List] := AllTrue[sc, MatchQ[{_, _, _, _?switchingCostValueQ}]];
switchingSpecQ[_] := False;

(* --- Validation --- *)

scenarioInputProblem[assoc_Association] :=
    Module[{model, graph, missing, vertices, adjacency, entries, exits, n,
            entryVertices, exitVertices, tolls, edgeSet},
        model = Lookup[assoc, "Model", Missing["KeyAbsent", "Model"]];
        If[MissingQ[model],
            Return[scenarioFailure["\"Model\" key is absent.", <|"MissingKeys" -> {"Model"}|>], Module]];
        If[!AssociationQ[model],
            Return[scenarioFailure["\"Model\" value must be an Association."], Module]];
        If[KeyExistsQ[assoc, "Data"],
            Return[scenarioFailure[
                "\"Data\" key is no longer supported. Provide numeric values directly in \"Model\"."], Module]];
        If[KeyExistsQ[model, "edgeEntryTolls"],
            Return[scenarioFailure["\"edgeEntryTolls\" was renamed to \"EdgeEntryTolls\"."], Module]];

        graph = Lookup[model, "Graph"];
        vertices = Lookup[model, "Vertices"];
        If[MatchQ[graph, _Graph] && ListQ[vertices] && !SubsetQ[VertexList[graph], vertices],
            Return[scenarioFailure["\"Vertices\" must all be vertices of \"Graph\"."], Module]];

        missing = Select[$RequiredModelKeys, !KeyExistsQ[model, #] &];
        If[KeyExistsQ[model, "EdgeEntryTolls"], missing = DeleteCases[missing, "Switching"]];
        If[missing =!= {},
            Return[scenarioFailure["Missing required Model keys: " <> StringRiffle[missing, ", "],
                <|"MissingKeys" -> missing|>], Module]];

        {adjacency, entries, exits} = Lookup[model, {"Adjacency", "Entries", "Exits"}];
        If[!(ListQ[vertices] && vertices =!= {} && AllTrue[vertices, IntegerQ] && DuplicateFreeQ[vertices]),
            Return[scenarioFailure["\"Vertices\" must be a nonempty list of distinct integers."], Module]];
        n = Length[vertices];
        If[!(MatrixQ[adjacency, NumericQ] && Dimensions[adjacency] === {n, n}),
            Return[scenarioFailure["\"Adjacency\" must be a numeric square matrix with one row per vertex."], Module]];
        If[!AllTrue[Normal[Diagonal[adjacency]], # == 0 &],
            Return[scenarioFailure["\"Adjacency\" must have a zero diagonal; self-loops are not supported."], Module]];

        If[!AllTrue[{entries, exits}, ListQ[#] && AllTrue[#, MatchQ[{_, _?realNumberQ}]] &],
            Return[scenarioFailure["\"Entries\" and \"Exits\" must be lists of {vertex, value} pairs with real values."], Module]];
        If[!AllTrue[entries[[All, 2]], # >= 0 &],
            Return[scenarioFailure["Entry inflows must be nonnegative."], Module]];
        entryVertices = First /@ entries;
        exitVertices  = First /@ exits;
        If[!SubsetQ[vertices, Join[entryVertices, exitVertices]],
            Return[scenarioFailure["Entry and exit vertices must belong to \"Vertices\"."], Module]];
        If[!DuplicateFreeQ[entryVertices] || !DuplicateFreeQ[exitVertices],
            Return[scenarioFailure["Each vertex may appear at most once in \"Entries\" and at most once in \"Exits\"."], Module]];
        If[Intersection[entryVertices, exitVertices] =!= {},
            Return[scenarioFailure["Entry and exit vertices must be disjoint."], Module]];

        If[!switchingSpecQ[Lookup[model, "Switching", {}]],
            Return[scenarioFailure["\"Switching\" must be {r, v, w, cost} rows or an Association {r, v, w} -> cost, with each cost a nonnegative number, Infinity, or a pure function."], Module]];

        edgeSet = modelEdgeSet[vertices, adjacency];
        tolls = Lookup[model, "EdgeEntryTolls", <||>];
        If[!(AssociationQ[tolls] && AllTrue[Keys[tolls], KeyExistsQ[edgeSet, #] &] &&
             AllTrue[Values[tolls], switchingCostValueQ]),
            Return[scenarioFailure["\"EdgeEntryTolls\" must be an Association {v, w} -> toll keyed by edges of the network, with each toll a nonnegative number, Infinity, or a pure function."], Module]];

        hamiltonianProblem[Lookup[assoc, "Hamiltonian", <||>], edgeSet]
    ];

hamiltonianProblem[spec_, edgeSet_Association] :=
    Module[{unknown, ham, edgeParams, badEdges, badValues},
        If[!AssociationQ[spec],
            Return[scenarioFailure["\"Hamiltonian\" must be an Association."], Module]];
        unknown = Complement[Keys[spec], Keys[$DefaultHamiltonian]];
        If[unknown =!= {},
            Return[scenarioFailure[
                "Unknown \"Hamiltonian\" keys: " <> StringRiffle[ToString /@ unknown, ", "] <>
                ". Allowed keys: " <> StringRiffle[Keys[$DefaultHamiltonian], ", "] <> ".",
                <|"UnknownKeys" -> unknown|>], Module]];
        ham = normalizeHamiltonianSpec[spec];
        If[!NumericQ[ham["Alpha"]],
            Return[scenarioFailure["\"Hamiltonian\" default \"Alpha\" must be numeric."], Module]];
        If[!hamiltonianTermQ[ham["V"]],
            Return[scenarioFailure["\"Hamiltonian\" default \"V\" must be numeric or a pure function."], Module]];
        If[!hamiltonianTermQ[ham["G"]],
            Return[scenarioFailure["\"Hamiltonian\" default \"G\" must be numeric or a pure function."], Module]];
        edgeParams = KeyTake[ham, {"EdgeAlpha", "EdgeV", "EdgeG"}];
        If[!AllTrue[edgeParams, AssociationQ],
            Return[scenarioFailure["\"Hamiltonian\" \"EdgeAlpha\", \"EdgeV\" and \"EdgeG\" must be Associations."], Module]];
        badEdges = Select[DeleteDuplicates[Join @@ (Keys /@ Values[edgeParams])], !KeyExistsQ[edgeSet, #] &];
        If[badEdges =!= {},
            Return[scenarioFailure["\"Hamiltonian\" edge-parameter keys must be edges {u,v} of the network.",
                <|"InvalidEdges" -> badEdges|>], Module]];
        badValues = Join[
            Select[Normal[ham["EdgeAlpha"]], !NumericQ[Last[#]] &],
            Select[Join[Normal[ham["EdgeV"]], Normal[ham["EdgeG"]]], !hamiltonianTermQ[Last[#]] &]];
        If[badValues =!= {},
            Return[scenarioFailure["\"Hamiltonian\" edge-parameter values must be numeric (or pure functions for EdgeV/EdgeG).",
                <|"InvalidRules" -> badValues|>], Module]];
        None
    ];

switchingProblem[model_Association, auxTriples_List] :=
    Module[{admissible, invalid},
        admissible = AssociationMap[True &, auxTriples];
        invalid = Select[Keys[normalizeSwitchingCosts[Lookup[model, "Switching", <||>]]],
            !KeyExistsQ[admissible, #] &];
        If[invalid === {},
            None,
            scenarioFailure["\"Switching\" keys must be admissible transitions {r, v, w}: {r, v} and {v, w} are edges and r != w.",
                <|"InvalidTriples" -> invalid|>]]
    ];

checkedTopology[assoc_Association] :=
    Module[{problem, model, topology},
        problem = scenarioInputProblem[assoc];
        If[FailureQ[problem], Return[problem, Module]];
        model = assoc["Model"];
        topology = Lookup[assoc, "Topology", buildAuxiliaryTopology[model]];
        If[!AssociationQ[topology],
            Return[scenarioFailure["Model topology is invalid or could not be constructed."], Module]];
        problem = switchingProblem[model, topology["AuxTriples"]];
        If[FailureQ[problem], problem, topology]
    ];

validateScenario::invalid = "Invalid scenario: `1`";
validateScenario::notscenario = "validateScenario expected a scenario object; got `1`.";

validateScenario[scenario[assoc_Association]] :=
    With[{topology = checkedTopology[assoc]},
        If[FailureQ[topology],
            Message[validateScenario::invalid, topology["Message"]]; False,
            True]
    ];
validateScenario[x_] := (Message[validateScenario::notscenario, Short[x]]; False);

(* --- Switching costs --- *)

expandEdgeEntryTolls[model_Association] :=
    Module[{tolls, vertices, neighbours},
        tolls = Lookup[model, "EdgeEntryTolls", <||>];
        If[tolls === <||>, Return[<||>, Module]];
        vertices   = model["Vertices"];
        neighbours = AssociationThread[vertices,
            Pick[vertices, #, 1] & /@ Normal[undirectedAdjacency[model["Adjacency"]]]];
        (* A toll on entering edge {v, w} at v applies to every transition
           {r, v, w} from another neighbour r of v; r = w would be the U-turn
           {w, v, w}. Numeric tolls become constant functions of the flow. *)
        Association @ Flatten @ KeyValueMap[
            Function[{edge, toll},
                With[{v = First[edge], w = Last[edge],
                      cost = If[MatchQ[toll, _Function], toll, Function[toll]]},
                    ({#, v, w} -> cost) & /@ DeleteCases[neighbours[v], w]]],
            tolls]
    ];

(* Exact Floyd-Warshall over the edge labels around one junction. GraphDistance
   on a weighted Graph coerces exact rational costs to machine reals, which
   would leak floats into the symbolic pipeline. *)
junctionClosure[sc_Association, triples_List] :=
    Module[{labels, index, n, dist},
        labels = DeleteDuplicates[Join[triples[[All, 1]], triples[[All, 3]]]];
        n      = Length[labels];
        index  = AssociationThread[labels, Range[n]];
        dist   = ConstantArray[Infinity, {n, n}];
        Do[dist[[k, k]] = 0, {k, n}];
        Scan[(dist[[index[#[[1]]], index[#[[3]]]]] = sc[#]) &, triples];
        Do[dist[[a, b]] = Min[dist[[a, b]], dist[[a, k]] + dist[[k, b]]], {k, n}, {a, n}, {b, n}];
        (# -> dist[[index[#[[1]]], index[#[[3]]]]]) & /@ triples
    ];

(* c(r, v, w) becomes the cheapest sequence of switches r -> k -> ... -> w at
   junction v. Infinity (blocked) and Function costs are neither updated nor
   used as intermediate steps. Keys keep their order. *)
metricClosureSwitchingCosts[sc_Association] :=
    Join[sc, Association @ Flatten @ Values @ GroupBy[
        Select[Keys[sc], NumericQ[sc[#]] &], #[[2]] &, junctionClosure[sc, #] &]];

(* Missing transitions cost 0 before the closure runs, so the finite costs
   the solver sees satisfy the triangle inequality
   c(r, v, w) <= c(r, v, k) + c(k, v, w) at every junction v. *)
completedSwitchingCosts[model_Association, auxTriples_List] :=
    metricClosureSwitchingCosts[Join[
        AssociationMap[0 &, auxTriples],
        expandEdgeEntryTolls[model],
        normalizeSwitchingCosts[Lookup[model, "Switching", <||>]]]];

(* --- Topology --- *)

auxTriplesFromPairs[directedPairs_List] :=
    Module[{incoming, outgoing},
        incoming = GroupBy[directedPairs, Last -> First];
        outgoing = GroupBy[directedPairs, First -> Last];
        DeleteDuplicates @ Flatten[
            Table[
                If[vIn =!= vOut, {vIn, vMid, vOut}, Nothing],
                {vMid, Intersection[Keys[incoming], Keys[outgoing]]},
                {vIn,  incoming[vMid]},
                {vOut, outgoing[vMid]}
            ],
            2
        ]
    ];

buildAuxiliaryTopology[model_Association] :=
    Module[{vertices, adjacency, entries, exits, graph, auxEntryVertices, auxExitVertices,
            entryEdges, exitEdges, halfPairs, inAuxEntryPairs, outAuxExitPairs, pairs, auxPairs},
        {vertices, adjacency, entries, exits} = Lookup[model, {"Vertices", "Adjacency", "Entries", "Exits"}];
        If[!ListQ[vertices] || !MatrixQ[adjacency] || !ListQ[entries] || !ListQ[exits],
            Return[$Failed, Module]];

        graph = Check[
            AdjacencyGraph[vertices, undirectedAdjacency[adjacency], DirectedEdges -> False],
            $Failed
        ];
        If[!MatchQ[graph, _Graph], Return[$Failed, Module]];

        auxEntryVertices = auxEntry /@ First /@ entries;
        auxExitVertices  = auxExit /@ First /@ exits;
        entryEdges = MapThread[DirectedEdge, {auxEntryVertices, First /@ entries}];
        exitEdges  = MapThread[DirectedEdge, {First /@ exits, auxExitVertices}];

        halfPairs       = List @@@ EdgeList[graph];
        inAuxEntryPairs = List @@@ entryEdges;
        outAuxExitPairs = List @@@ exitEdges;
        pairs           = Join[halfPairs, Reverse[halfPairs, {2}]];
        auxPairs        = Join[inAuxEntryPairs, outAuxExitPairs, pairs];

        <|
            "Graph"             -> graph,
            "AuxiliaryGraph"    -> EdgeAdd[graph, Join[entryEdges, exitEdges]],
            "AuxEntryVertices"  -> auxEntryVertices,
            "AuxExitVertices"   -> auxExitVertices,
            "AuxEntryEdges"     -> entryEdges,
            "AuxExitEdges"      -> exitEdges,
            "AuxTriples"        -> auxTriplesFromPairs[auxPairs],
            "HalfPairs"         -> halfPairs,
            "InAuxEntryPairs"   -> inAuxEntryPairs,
            "OutAuxExitPairs"   -> outAuxExitPairs,
            "Pairs"             -> pairs,
            "AuxPairs"          -> auxPairs
        |>
    ];

(* --- Type predicate and accessor --- *)

scenarioQ[x_] := mfgTypedQ[x, scenario];

scenarioData[s_] := mfgData[s];
scenarioData[s_, key_] := mfgData[s, key];

(* --- Complete --- *)

completeScenario::notscenario =
    "completeScenario expected a scenario object; got `1`. Returning input unchanged.";
completeScenario::notopology =
    "completeScenario was called without cached topology; rebuilding it from the Model and storing it in the scenario.";

completeScenario[scenario[assoc_Association]] :=
    Module[{model, topology},
        model    = Lookup[assoc, "Model", <||>];
        topology = Lookup[assoc, "Topology", Missing["KeyAbsent", "Topology"]];
        If[MissingQ[topology],
            Message[completeScenario::notopology];
            topology = If[AssociationQ[model], buildAuxiliaryTopology[model], $Failed]
        ];
        scenario[Join[
            assoc,
            <|
                "Hamiltonian" -> normalizeHamiltonianSpec[Lookup[assoc, "Hamiltonian", <||>]],
                "Identity"    -> Lookup[assoc, "Identity", <||>],
                "Benchmark"   -> Join[
                    <|"Tier" -> $DefaultBenchmarkTier, "Timeout" -> $DefaultBenchmarkTimeout|>,
                    Lookup[assoc, "Benchmark", <||>]]
            |>,
            If[AssociationQ[topology],
                <|"Topology" -> topology,
                  "Model"    -> Join[model, <|"Switching" -> completedSwitchingCosts[model, topology["AuxTriples"]]|>]|>,
                <||>
            ]
        ]]
    ];

completeScenario[x_] := (Message[completeScenario::notscenario, x]; x);

(* --- Constructor --- *)

makeScenario[rawAssoc_Association] :=
    Module[{rawModel, provenance, assoc, topology},
        rawModel = Lookup[rawAssoc, "Model", Missing["KeyAbsent", "Model"]];

        (* Preserve mathematically relevant source data before normalization,
           metric closure, completion, or other derived transformations can
           erase evidence of inexact user input. *)
        provenance = If[AssociationQ[rawModel],
            <|
                "Model" -> KeyTake[rawModel, {"Entries", "Exits", "Switching", "EdgeEntryTolls"}],
                "Hamiltonian" -> Lookup[rawAssoc, "Hamiltonian", <||>]
            |>,
            <||>
        ];
        assoc = If[AssociationQ[rawModel],
            Join[rawAssoc, <|"Model" -> normalizeScenarioModel[rawModel]|>],
            rawAssoc
        ];

        (* Always rebuild the topology from the model, never reuse one passed in. *)
        topology = checkedTopology[KeyDrop[assoc, "Topology"]];
        If[FailureQ[topology], Return[topology, Module]];

        completeScenario[scenario[Join[assoc,
            <|"Topology" -> topology, "InputProvenance" -> provenance|>]]]
    ];

makeScenario[_] :=
    scenarioFailure["makeScenario requires an Association as input."];

End[]

Protect[auxEntry, auxExit];

EndPackage[]
