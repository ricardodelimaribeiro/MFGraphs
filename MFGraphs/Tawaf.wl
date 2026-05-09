(* Wolfram Language package *)
(* Tawaf: specialized scenario builder for Tawaf circumambulation with shared
   physical-edge congestion across rounds.

   Modelling note: docs/research/notes/tawaf-model.md
   - ritual context, "unroll then couple" idea, scope of coupling
     (which system blocks are rewritten and why), parameter provenance
     (Black Stone potential, entry flow, layer interpretation). *)

BeginPackage["Tawaf`",
    {"primitives`", "utilities`", "scenarioTools`",
     "unknownsTools`", "systemTools`"}];

makeTawafScenario::usage =
"makeTawafScenario[rounds, nodesPerRound] and makeTawafScenario[rounds, nodesPerRound, layers] build an unrolled Tawaf scenario with rounds-by-nodesPerRound logical nodes per layer (default 1 layer). Tangential edges go forward only (counter-clockwise); radial edges (when layers > 1) connect adjacent layers in both directions. Black-Stone potential V = -5 is set on edges leaving position 1 of any round. Construction parameters are stored in scenarioData[s, \"Tawaf\"].";

makeTawafSystem::usage =
"makeTawafSystem[s] and makeTawafSystem[s, unk] build an mfgSystem and rewrite EqGeneral and AltOptCond so that flows on the same physical edge (same position-pair-and-direction) share congestion. Construction dimensions are read from scenarioData[s, \"Tawaf\"]; returns a Failure if that metadata is missing.";

Begin["`Private`"];

(* ============================================================
   Coordinate codec: single source of truth.
   Layer-major flattening: vertex = (l-1)*rounds*nodesPerRound + (r-1)*nodesPerRound + p
   ============================================================ *)

tawafEncode[r_Integer, p_Integer, l_Integer, rounds_Integer, nodesPerRound_Integer] :=
    (l - 1)*rounds*nodesPerRound + (r - 1)*nodesPerRound + p;

(* Returns {round, position, layer} (1-based). *)
tawafDecode[v_Integer, rounds_Integer, nodesPerRound_Integer] :=
    Module[{layerStride = rounds*nodesPerRound, idx, l, rem, r, p},
        idx = v - 1;
        l = Quotient[idx, layerStride] + 1;
        rem = Mod[idx, layerStride];
        r = Quotient[rem, nodesPerRound] + 1;
        p = Mod[rem, nodesPerRound] + 1;
        {r, p, l}
    ];

(* Physical-edge identifier for a logical flow variable j[u, v].
   - Both endpoints must be integers (logical nodes). Any non-integer endpoint
     (e.g. "auxEntry1" / "auxExit3" strings) marks the edge as boundary.
   - Tangential: same layer; key by ordered position pair, layer, direction.
   - Radial: same position; key by position, layer pair, direction. *)
tawafPhysicalId[var_, rounds_Integer, nodesPerRound_Integer] :=
    With[{a = var[[1]], b = var[[2]]},
        Which[
            !IntegerQ[a] || !IntegerQ[b], "Boundary",
            True,
                Module[{ra, pa, la, rb, pb, lb},
                    {ra, pa, la} = tawafDecode[a, rounds, nodesPerRound];
                    {rb, pb, lb} = tawafDecode[b, rounds, nodesPerRound];
                    Which[
                        la === lb,
                            (* Tangential motion within a layer.
                               Direction: forward iff pb == pa mod nodesPerRound + 1
                               (handles wrap-around pa = nodesPerRound, pb = 1). *)
                            With[{isForward = (pb == Mod[pa, nodesPerRound] + 1)},
                                {"Tangential", Min[pa, pb], Max[pa, pb], la, isForward}
                            ],
                        pa === pb,
                            (* Radial motion between layers at the same position. *)
                            With[{isOutward = lb > la},
                                {"Radial", pa, Min[la, lb], Max[la, lb], isOutward}
                            ],
                        True, "Boundary"  (* unrecognised topology *)
                    ]
                ]
        ]
    ];

(* Build coupling rules: each logical flow is replaced by the sum of all
   logical flows mapping to the same physical id. Singleton groups and the
   "Boundary" group are skipped (no rewrite needed). *)
tawafCouplingRules[js_List, rounds_Integer, nodesPerRound_Integer] :=
    Module[{groups},
        groups = GroupBy[js, tawafPhysicalId[#, rounds, nodesPerRound] &];
        Flatten @ KeyValueMap[
            Function[{id, vars},
                If[id === "Boundary" || Length[vars] < 2,
                    {},
                    With[{total = Total[vars]}, (# -> total) & /@ vars]
                ]
            ],
            groups
        ]
    ];

(* Helper: rebuild a typed sub-record by replacing one key. *)
replaceInTyped[record_, head_, key_String, newValue_] :=
    head[Join[mfgData[record], <|key -> newValue|>]];

(* ============================================================
   Scenario builder
   ============================================================ *)

makeTawafScenario[rounds_Integer, nodesPerRound_Integer] :=
    makeTawafScenario[rounds, nodesPerRound, 1];

makeTawafScenario[rounds_Integer, nodesPerRound_Integer, layers_Integer] /;
    rounds >= 1 && nodesPerRound >= 2 && layers >= 1 :=
    Module[{enc, totalNodes, tangentialEdges, radialEdges, allEdges, graph,
            entries, exits, blackStoneNodes, vMap, edgeV, raw, scenarioObj},

        enc[r_, p_, l_] := tawafEncode[r, p, l, rounds, nodesPerRound];
        totalNodes = rounds*nodesPerRound*layers;

        (* Forward-only tangential edges: counter-clockwise circumambulation. *)
        tangentialEdges = Flatten @ Table[
            If[p < nodesPerRound,
                DirectedEdge[enc[r, p, l], enc[r, p + 1, l]],
                If[r < rounds,
                    DirectedEdge[enc[r, p, l], enc[r + 1, 1, l]],
                    Nothing
                ]
            ],
            {l, layers}, {r, rounds}, {p, nodesPerRound}
        ];

        (* Radial edges: bidirectional between adjacent layers. *)
        radialEdges = Flatten @ Table[
            {DirectedEdge[enc[r, p, l],     enc[r, p, l + 1]],
             DirectedEdge[enc[r, p, l + 1], enc[r, p, l]]},
            {r, rounds}, {p, nodesPerRound}, {l, layers - 1}
        ];

        allEdges = Join[tangentialEdges, radialEdges];
        graph = Graph[Range[totalNodes], allEdges, DirectedEdges -> True];

        (* Boundary: entries at first node of round 1 (each layer); exits at
           the last node of the last round (each layer). *)
        entries = Table[{enc[1, 1, l], 100/layers}, {l, layers}];
        exits   = Table[{enc[rounds, nodesPerRound, l], 0}, {l, layers}];

        (* Black-Stone potential: -5 on edges originating at position 1. *)
        blackStoneNodes = Flatten @ Table[enc[r, 1, l], {r, rounds}, {l, layers}];
        vMap = AssociationThread[blackStoneNodes, -5.0];
        edgeV = Association @ Map[
            (List @@ #) -> Lookup[vMap, First[#], 0.0] &,
            allEdges
        ];

        raw = <|
            "Model" -> <|
                "Graph"     -> graph,
                "Entries"   -> entries,
                "Exits"     -> exits,
                "Switching" -> {}
            |>,
            "Hamiltonian" -> <|"EdgeV" -> edgeV|>,
            "Tawaf" -> <|
                "Rounds"        -> rounds,
                "NodesPerRound" -> nodesPerRound,
                "Layers"        -> layers
            |>
        |>;

        scenarioObj = makeScenario[raw];
        If[FailureQ[scenarioObj], Return[scenarioObj, Module]];
        scenarioObj
    ];

makeTawafScenario[rounds_, nodesPerRound_, layers_] :=
    Failure["makeTawafScenario", <|
        "MessageTemplate"   -> "Invalid Tawaf dimensions: rounds=`1`, nodesPerRound=`2`, layers=`3` (must be positive integers; nodesPerRound >= 2).",
        "MessageParameters" -> {rounds, nodesPerRound, layers}
    |>];

(* ============================================================
   System builder
   ============================================================ *)

makeTawafSystem[s_?scenarioQ] :=
    Module[{unk = makeSymbolicUnknowns[s]},
        makeTawafSystem[s, unk]
    ];

makeTawafSystem[s_?scenarioQ, unk_?symbolicUnknownsQ] :=
    Module[{meta, rounds, nodesPerRound, layers, sys, js, rules,
            data, ham, comp},
        meta = scenarioData[s, "Tawaf"];
        If[!AssociationQ[meta] ||
           !AllTrue[{"Rounds", "NodesPerRound", "Layers"}, KeyExistsQ[meta, #] &],
            Return[Failure["makeTawafSystem", <|
                "MessageTemplate"   -> "Scenario does not carry Tawaf metadata; cannot derive (rounds, nodesPerRound, layers).",
                "MessageParameters" -> {}
            |>], Module]
        ];
        rounds        = meta["Rounds"];
        nodesPerRound = meta["NodesPerRound"];
        layers        = meta["Layers"];

        sys = makeSystem[s, unk];
        js  = systemData[sys, "Js"];
        rules = tawafCouplingRules[js, rounds, nodesPerRound];

        If[rules === {}, Return[sys, Module]];

        data = mfgData[sys];
        ham  = data["HamiltonianData"];
        comp = data["ComplementarityData"];

        ham = replaceInTyped[ham, mfgHamiltonianData, "EqGeneral",
            mfgData[ham]["EqGeneral"] /. rules];
        comp = replaceInTyped[comp, mfgComplementarityData, "AltOptCond",
            mfgData[comp]["AltOptCond"] /. rules];

        mfgSystem[Join[data, <|
            "HamiltonianData"      -> ham,
            "ComplementarityData"  -> comp
        |>]]
    ];

End[];

EndPackage[];
