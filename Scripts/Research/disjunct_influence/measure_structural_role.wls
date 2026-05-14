#!/usr/bin/env wolframscript
(* PR 2: per-abstract-complementarity structural features + outcome labels.

   For each scenario:
     1. Enumerate abstract complementarities from systemData's Alt* buckets
        (AltFlows, AltTransitionFlows, AltOptCond, AltExitCond). One row
        per Or expression.
     2. Compute per-disjunct features:
          - anchorVar      : the j-variable inside the Or (first _j match)
          - anchorVertex   : middle index (for j[a,v,b]) or head (for j[a,b])
          - anchorEdge     : (a,v) or (a,b) underlying directed edge
          - junctionDegree : in+out degree of anchorVertex in the auxiliary graph
          - boundaryDistance: BFS distance from anchorVertex to nearest entry/exit
          - edgeBetweenness: per-edge betweenness centrality from
                             EdgeBetweennessCentrality (Wolfram built-in)
          - edgeCost       : Lookup[scenario["Switching"], anchorEdge, 0]
     3. Outcome label from solutions/<key>.wxf:
          - Inactive : rule var -> 0 (or |val| < 10^-8)
          - Active   : rule var -> v with |v| > 10^-8
          - Ambiguous: var appears in Residual, OR not present in Rules,
                       OR no cached solution
     4. Emit JSON to Scripts/Research/results/structural/<key>.json. *)

$RepoRoot = ParentDirectory[ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]]];
PrependTo[$Path, FileNameJoin[{$RepoRoot, "MFGraphs"}]];
Needs["MFGraphs`"];
Get[FileNameJoin[{$RepoRoot, "Scripts", "SolutionCacheHelpers.wls"}]];

$resultsDir = FileNameJoin[{$RepoRoot, "Scripts", "Research", "results", "structural"}];
If[!DirectoryQ[$resultsDir], CreateDirectory[$resultsDir, CreateIntermediateDirectories -> True]];

sanitizeKey[k_Integer] := "case_" <> ToString[k];
sanitizeKey[k_String]  := StringReplace[k, " " -> "_"];

(* Anchor extraction: grab the first j-variable from the Or; its argument
   tuple identifies the anchor vertex and edge. *)
extractAnchor[orExpr_] :=
    Module[{js},
        js = Cases[orExpr, _j, {0, Infinity}];
        If[js === {}, Return[<|"var" -> Missing[], "vertex" -> Missing[], "edge" -> Missing[]|>]];
        With[{firstJ = First[js], args = List @@ First[js]},
            Switch[Length[args],
                2, <|"var" -> firstJ, "vertex" -> args[[2]], "edge" -> args|>,
                3, <|"var" -> firstJ, "vertex" -> args[[2]], "edge" -> {args[[1]], args[[2]]}|>,
                _, <|"var" -> firstJ, "vertex" -> Missing[], "edge" -> Missing[]|>
            ]
        ]
    ];

(* Outcome label from cached solution. *)
$tightThreshold = 10^-8;

outcomeLabel[var_, cached_Association] :=
    Module[{sol, ruleSet, residual, value},
        sol = Lookup[cached, "Solution", Missing[]];
        (* Solution shapes seen in solutions/*.wxf:
             - <|"Rules" -> {var -> v, ...}, "Residual" -> True | parametric Or|>
             - flat <|var -> v, var -> v, ...|>     (fully-determined case)
             - flat list of Rule[var, v]            (older cache format) *)
        ruleSet = Which[
            AssociationQ[sol] && KeyExistsQ[sol, "Rules"], sol["Rules"],
            AssociationQ[sol],                              sol,
            ListQ[sol] && AllTrue[sol, MatchQ[#, _Rule | _RuleDelayed] &], Association[sol],
            True,                                           Missing[]
        ];
        residual = If[AssociationQ[sol], Lookup[sol, "Residual", Null], Null];
        If[residual =!= Null && residual =!= True && !FreeQ[residual, var],
            Return["Ambiguous"]];
        If[ruleSet === Missing[], Return["Ambiguous"]];
        ruleSet = If[ListQ[ruleSet], Association[ruleSet], ruleSet];
        If[!KeyExistsQ[ruleSet, var], Return["Ambiguous"]];
        value = ruleSet[var];
        Which[
            !NumericQ[value],             "Ambiguous",
            Abs[value] < $tightThreshold, "Inactive",
            True,                         "Active"
        ]
    ];
outcomeLabel[_, _] := "Ambiguous";

(* RawJSON serializer can't handle Missing[], Infinity, or symbolic heads.
   Walk the rows and replace those with Null / strings. *)
sanitizeForJSON[expr_] := expr //. {
    _Missing -> Null,
    Infinity -> Null,
    -Infinity -> Null,
    s_Symbol /; !MemberQ[{True, False, Null}, s] :> ToString[s, InputForm]
};

measureScenario[label_String, key_, scenario_] :=
    Module[{sys, topo, graph, entryVerts, exitVerts, allBoundary, edgeList,
            altOrs, anchors, betweenness, betweennessAssoc, vertexDegrees,
            switching, cached, rows, json, outFile, jsonRows},
        sys = makeSystem[scenario];
        topo = scenarioData[scenario, "Topology"];
        graph = topo["AuxiliaryGraph"];
        entryVerts = First /@ scenarioData[scenario, "Model"]["Entries"];
        exitVerts  = First /@ scenarioData[scenario, "Model"]["Exits"];
        allBoundary = Join[entryVerts, exitVerts];
        (* VertexDegree on graphs whose VertexList contains a mix of
           Integer and String vertices crashes the kernel ("license error").
           Compute degree from EdgeList ourselves to dodge the bug. *)
        edgeList = EdgeList[graph];
        vertexDegrees = Counts[Flatten[List @@@ edgeList]];
        switching = Replace[
            systemData[sys, "SwitchingCosts"],
            Except[_Association] -> <||>
        ];
        cached = Quiet @ loadSolution[key];
        If[!AssociationQ[cached], cached = <||>];
        (* Collect all Or expressions from Alt* buckets. *)
        altOrs = Cases[
            Flatten[{
                systemData[sys, "AltFlows"],
                systemData[sys, "AltTransitionFlows"],
                systemData[sys, "AltOptCond"],
                systemData[sys, "AltExitCond"]
            }],
            _Or
        ];
        (* Edge betweenness centrality, Association keyed by directed edge. *)
        (* Build betweenness lookup keyed by an unordered {head,tail} list
           to dodge head-display ambiguities (DirectedEdge vs UndirectedEdge
           Lookup mismatches we observed on graphs with mixed-type vertices). *)
        betweennessAssoc = Module[{vals = Quiet @ Check[EdgeBetweennessCentrality[graph], {}]},
            If[Length[vals] === Length[edgeList],
                Association @ MapThread[
                    Function[{e, v},
                        Sort[List @@ e] -> v
                    ],
                    {edgeList, vals}],
                <||>]
        ];
        anchors = extractAnchor /@ altOrs;
        rows = MapThread[
            Function[{orExpr, anchor},
                Module[{vertex = anchor["vertex"], edge = anchor["edge"], var = anchor["var"],
                        outcome, junctionDegree, boundaryDistance, edgeBetweenness, edgeCost,
                        directedEdge},
                    junctionDegree = Lookup[vertexDegrees, vertex, Missing[]];
                    boundaryDistance = If[(IntegerQ[vertex] || StringQ[vertex]) && allBoundary =!= {},
                        Module[{ds = Quiet @ Check[GraphDistance[graph, vertex, #], Infinity] & /@ allBoundary},
                            With[{finite = Cases[ds, _?(NumericQ[#] && # < Infinity &)]},
                                If[finite === {}, Missing[], Min[finite]]
                            ]
                        ],
                        Missing[]];
                    edgeBetweenness = If[ListQ[edge] && Length[edge] === 2,
                        Lookup[betweennessAssoc, Key[Sort[edge]], Missing[]],
                        Missing[]];
                    (* SwitchingCosts is keyed by 3-tuples {r,i,w}; edge here may
                       be a 2- or 3-tuple. Look up directly only for matching
                       arity; otherwise sum over all 3-tuples that include this
                       2-edge as either incoming or outgoing leg. *)
                    edgeCost = If[ListQ[edge],
                        Switch[Length[edge],
                            3, Lookup[switching, Key[edge], 0],
                            2, Total[Values[KeySelect[switching,
                                MatchQ[#, {Sequence @@ edge, _} | {_, Sequence @@ edge}] &]]],
                            _, Missing[]
                        ],
                        Missing[]];
                    outcome = If[var === Missing[], "Ambiguous", outcomeLabel[var, cached]];
                    <|
                        "anchorVar"        -> ToString[var, InputForm],
                        "anchorVertex"    -> ToString[vertex],
                        "anchorEdge"      -> ToString[edge, InputForm],
                        "junctionDegree" -> If[NumericQ[junctionDegree], junctionDegree, Null],
                        "boundaryDistance" -> If[NumericQ[boundaryDistance] && boundaryDistance < Infinity, boundaryDistance, Null],
                        "edgeBetweenness" -> If[NumericQ[edgeBetweenness], N[edgeBetweenness], Null],
                        "edgeCost"        -> If[NumericQ[edgeCost], N[edgeCost], Null],
                        "outcome"         -> outcome
                    |>
                ]
            ],
            {altOrs, anchors}
        ];
        json = <|
            "scenario"     -> label,
            "key"          -> ToString[key],
            "abstractCount" -> Length[rows],
            "hasCachedSolution" -> AssociationQ[cached] && cached =!= <||>,
            "outcomeCounts" -> Association[
                Function[lbl, lbl -> Count[rows, KeyValuePattern["outcome" -> lbl]]] /@
                {"Inactive", "Active", "Ambiguous"}],
            "rows"          -> rows
        |>;
        outFile = FileNameJoin[{$resultsDir, label <> ".json"}];
        Export[outFile, json, "JSON"];
        Print[
            StringPadRight[label, 20],
            " abstract=", Length[rows],
            " cached=", AssociationQ[cached] && cached =!= <||>,
            " outcomes=", Counts[#["outcome"] & /@ rows]
        ]
    ];

Print["Sweep: ", Length[$ExampleBoundaryDefaults], " scenarios, results to ", $resultsDir];
Scan[
    Function[k,
        Module[{ee = $ExampleBoundaryDefaults[k]},
            Print["[start] ", k];
            Quiet @ Check[
                measureScenario[sanitizeKey[k], k, getExampleScenario[k, ee[[1]], ee[[2]]]],
                Print["[FAILED] ", k]
            ];
            Print["[done]  ", k]
        ]
    ],
    Keys[$ExampleBoundaryDefaults]
];
Print["Done."];
