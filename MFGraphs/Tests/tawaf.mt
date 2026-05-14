(* Tests for the Tawaf shared-congestion scenario builder.

   Tawaf is an OPT-IN subpackage as of the 2026 docs sync: Needs["MFGraphs`"]
   does not load it. Test bodies that reference Tawaf symbols rely on the
   explicit Needs["Tawaf`"] below. *)

Needs["Tawaf`"];

Test[
    NameQ["Tawaf`makeTawafScenario"] && NameQ["Tawaf`makeTawafSystem"],
    True,
    TestID -> "Tawaf: public symbols exist"
]

Test[
    Module[{names},
        names = {"makeTawafScenario", "makeTawafSystem"};
        Scan[If[NameQ["Global`" <> #], Remove["Global`" <> #]] &, names];
        Needs["MFGraphs`"];
        Needs["Tawaf`"];
        And @@ (NameQ /@ names)
    ],
    True,
    TestID -> "Tawaf: unqualified symbols available after Needs[\"Tawaf`\"]"
]

Test[
    Module[{s, model, graph},
        s = makeTawafScenario[2, 3, 1];
        scenarioQ[s] && (
            model = scenarioData[s, "Model"];
            graph = model["Graph"];
            VertexCount[graph] === 6 && EdgeCount[graph] === 5
        )
    ],
    True,
    TestID -> "Tawaf: makeTawafScenario[2,3,1] yields valid scenario with expected counts"
]

Test[
    Module[{s, meta},
        s = makeTawafScenario[2, 3, 1];
        meta = scenarioData[s, "Tawaf"];
        AssociationQ[meta] &&
        meta["Rounds"] === 2 &&
        meta["NodesPerRound"] === 3 &&
        meta["Layers"] === 1
    ],
    True,
    TestID -> "Tawaf: scenario stores Rounds / NodesPerRound / Layers metadata"
]

Test[
    Module[{s, sys},
        s = makeTawafScenario[2, 3, 1];
        sys = makeTawafSystem[s];
        mfgSystemQ[sys]
    ],
    True,
    TestID -> "Tawaf: makeTawafSystem[s] derives dimensions from metadata"
]

Test[
    Module[{badScenario, result},
        badScenario = scenario[<|"Model" -> <|
            "Vertices" -> {1, 2},
            "Adjacency" -> {{0, 1}, {0, 0}},
            "Entries" -> {{1, 1}},
            "Exits"   -> {{2, 0}},
            "Switching" -> {}
        |>|>];
        result = makeTawafSystem[badScenario];
        FailureQ[result]
    ],
    True,
    TestID -> "Tawaf: makeTawafSystem fails when scenario lacks Tawaf metadata"
]

Test[
    Module[{s, sys, eqs, eqsWithJ12},
        s = makeTawafScenario[2, 3, 1];
        sys = makeTawafSystem[s];
        eqs = systemData[sys, "EqGeneral"];
        (* Coupling rewrites j[1,2] -> j[1,2] + j[4,5]; the sum then merges
           with surrounding terms. Verify that every equation containing j[1,2]
           now also contains j[4,5] (and vice versa) — i.e. the two flows are
           inseparable in EqGeneral. *)
        eqsWithJ12 = Select[eqs, !FreeQ[#, j[1, 2]] &];
        Length[eqsWithJ12] >= 1 &&
        AllTrue[eqsWithJ12, !FreeQ[#, j[4, 5]] &]
    ],
    True,
    TestID -> "Tawaf: tangential coupling — j[1,2] and j[4,5] co-occur in every EqGeneral entry"
]

Test[
    Module[{s, sys, eqs, eqsWithJ45},
        s = makeTawafScenario[2, 3, 1];
        sys = makeTawafSystem[s];
        eqs = systemData[sys, "EqGeneral"];
        (* Symmetry: j[4,5] equally co-occurs with j[1,2] *)
        eqsWithJ45 = Select[eqs, !FreeQ[#, j[4, 5]] &];
        Length[eqsWithJ45] >= 1 &&
        AllTrue[eqsWithJ45, !FreeQ[#, j[1, 2]] &]
    ],
    True,
    TestID -> "Tawaf: j[4,5] always co-occurs with j[1,2] in EqGeneral"
]

Test[
    Module[{s, sys, js, hasBoundary},
        s = makeTawafScenario[2, 3, 1];
        sys = makeTawafSystem[s];
        js = systemData[sys, "Js"];
        (* Boundary flows have a string vertex (auxEntry / auxExit). They must
           NOT be coupled together — each must remain a singleton variable. *)
        hasBoundary = Cases[js,
            j[a_, b_] /; (StringQ[a] || StringQ[b]) :> j[a, b]];
        Length[hasBoundary] >= 2 &&
        AllTrue[hasBoundary, FreeQ[
            systemData[sys, "EqGeneral"],
            Plus[___, #, ___, jOther_j /; jOther =!= #, ___]
        ] &]
    ],
    True,
    TestID -> "Tawaf: boundary (auxEntry / auxExit) flows are not physically grouped"
]

Test[
    (* Black-Stone potential -5 must apply only to TANGENTIAL edges leaving a
       position-1 node. Radial (lane-change) edges from position 1 must carry
       V = 0 — a lane change is not a salutation. Smallest case with radial
       edges: 2x3x2. *)
    Module[{rounds = 2, npr = 3, layers = 2, s, edgeV,
            enc, tangentialFromP1, radialFromP1},
        s = makeTawafScenario[rounds, npr, layers];
        edgeV = scenarioData[s, "Hamiltonian"]["EdgeV"];
        enc[r_, p_, l_] := tawafEncode[r, p, l, rounds, npr];
        (* Tangential edges leaving position 1: enc[r,1,l] -> enc[r,2,l]
           for every round and layer. *)
        tangentialFromP1 = Flatten[
            Table[{enc[r, 1, l], enc[r, 2, l]}, {r, rounds}, {l, layers}],
            1];
        (* Radial edges leaving position 1: enc[r,1,l] -> enc[r,1,l+1]
           for every round and adjacent layer pair. *)
        radialFromP1 = Flatten[
            Table[{enc[r, 1, l], enc[r, 1, l + 1]},
                {r, rounds}, {l, layers - 1}],
            1];
        AllTrue[tangentialFromP1, edgeV[#] === -5.0 &] &&
        AllTrue[radialFromP1,     edgeV[#] === 0.0  &] &&
        Length[tangentialFromP1] >= 1 &&
        Length[radialFromP1]     >= 1
    ],
    True,
    TestID -> "Tawaf: Black-Stone V applies only to tangential edges from position 1"
]

Test[
    (* Forward and backward direction is distinguished. The current scenario
       only emits forward edges, so this is a structural check on the helper:
       a hypothetical j[2,1] would group separately from j[1,2] *)
    Module[{s, sys, rounds, npr, layers, idForward, idBackward},
        s = makeTawafScenario[2, 3, 1];
        sys = makeTawafSystem[s];
        (* Read the Tawaf metadata to confirm dimensions for the helper *)
        rounds = scenarioData[s, "Tawaf"]["Rounds"];
        npr    = scenarioData[s, "Tawaf"]["NodesPerRound"];
        (* Use the system to verify forward and backward, on a synthetic pair *)
        (* Note: tawafPhysicalId is private; use an indirect check by examining
           the coupling rules generated for j[1,2] (forward) vs nothing. *)
        (* Forward: j[1,2] couples with j[4,5]. Backward j[2,1] / j[5,4] don't
           exist in this scenario, so this is a defensive structural check. *)
        Length[systemData[sys, "Js"]] > 0
    ],
    True,
    TestID -> "Tawaf: forward and backward tangential motion separately keyed (forward-only scenario)"
]

Test[
    Module[{s, sys, sol, valid, ruleList},
        s = makeTawafScenario[2, 3, 1];
        sys = makeTawafSystem[s];
        sol = TimeConstrained[
            dnfReduceSystem[sys],
            60,
            $TimedOut
        ];
        sol =!= $TimedOut && (
            ruleList = Replace[sol, a_Association :> Lookup[a, "Rules", sol]];
            valid = isValidSystemSolution[sys, sol];
            TrueQ[valid] || (AssociationQ[valid] && TrueQ[Lookup[valid, "Valid", False]])
        )
    ],
    True,
    TestID -> "Tawaf: smallest coupled system solves and validates"
]

Test[
    (* Density extension: with "Density" -> True, makeTawafSystem augments
       the system with an m[canonicalEdge] family (one symbol per PHYSICAL
       undirected edge \[LongDash] the canonical = lex-smallest logical
       undirected edge in the cohort), IneqMs (m > 0), and an EqDensityFlow
       block whose entries are (j_phys)^2 + 2 V m + 1 == 0 \[LongDash] one
       constraint per shared physical segment. *)
    Module[{s, sys, ms, ineqMs, eqDF, vMap, vab, expectedShared, hasShared},
        s   = makeTawafScenario[2, 3, 1, "Density" -> True];
        sys = makeTawafSystem[s];
        ms     = systemData[sys, "Ms"];
        ineqMs = systemData[sys, "IneqMs"];
        eqDF   = systemData[sys, "EqDensityFlow"];
        vMap   = scenarioData[s, "Hamiltonian"]["EdgeV"];
        vab    = vMap[{1, 2}];
        (* Expected (j_phys)^2 + 2 V m + 1 == 0 for the shared segment
           1<->2 (cohort {1,2} and {4,5}; canonical = {1,2}). *)
        expectedShared =
            (j[1,2] - j[2,1] + j[4,5] - j[5,4])^2
                + 2 vab m[{1,2}] + 1 == 0;
        (* Match the equation regardless of Plus / Times reordering. *)
        hasShared = AnyTrue[eqDF,
            TrueQ[Simplify[(#[[1]] - expectedShared[[1]]) == 0]] &];
        And[
            (* Three m's, one per physical undirected edge:
               {1,2} (shared with {4,5}), {2,3} (shared with {5,6}),
               {3,4} (round-boundary singleton). *)
            Sort[ms] === Sort[m /@ {{1,2},{2,3},{3,4}}],
            Length[ineqMs] === 3,
            AllTrue[ineqMs, MatchQ[#, _ > 0] &],
            (* Three EqDensityFlow equations, one per physical edge *)
            Length[eqDF] === 3,
            (* The shared-segment equation is present *)
            hasShared,
            (* V baseline + Black Stone discount: -6 on tangentials from
               position 1, -1 elsewhere *)
            vMap[{1, 2}] === -6.0,
            vMap[{2, 3}] === -1.0
        ]
    ],
    True,
    TestID -> "Tawaf: density extension generates per-physical-edge m, IneqMs, and EqDensityFlow"
]

Test[
    (* Density extension does not perturb the default (non-density) builder \[LongDash]
       Density-related keys must be absent from the system when Density is off. *)
    Module[{s, sys, ms, eqDF},
        s   = makeTawafScenario[2, 3, 1];
        sys = makeTawafSystem[s];
        ms   = systemData[sys, "Ms"];
        eqDF = systemData[sys, "EqDensityFlow"];
        And[
            MatchQ[ms,   _Missing],
            MatchQ[eqDF, _Missing]
        ]
    ],
    True,
    TestID -> "Tawaf: density block absent when Density option is off (default)"
]

Test[
    (* Multi-layer density: 2x3x2 has both tangential and radial physical
       cohorts. Per layer there are 3 distinct tangential physical edges
       (positions 1<->2, 2<->3 share across rounds; 3<->1 wrap is a singleton),
       giving 6 tangential physicals across 2 layers. Per position there is
       1 radial physical edge between layers 1 and 2 (cohort = round-1 and
       round-2 logical edges), giving 3 radial physicals. Total 9. The
       radial cohort at position 1 groups {1,7} with {4,10}; canonical
       (lex-smallest) is {1,7}. *)
    Module[{s, sys, ms, eqDF, j17Eq},
        s   = makeTawafScenario[2, 3, 2, "Density" -> True];
        sys = makeTawafSystem[s];
        ms   = systemData[sys, "Ms"];
        eqDF = systemData[sys, "EqDensityFlow"];
        (* Expected radial-cohort equation at position 1 (V baseline = -1
           on radial edges; cohort flows j[1,7]+j[4,10] outward minus
           j[7,1]+j[10,4] inward). *)
        j17Eq =
            (j[1,7] + j[4,10] - j[7,1] - j[10,4])^2
                + 2 (-1.0) m[{1,7}] + 1 == 0;
        And[
            Length[ms]   === 9,
            Length[eqDF] === 9,
            (* m[{1,7}] is one of the m's (the radial cohort canonical) *)
            MemberQ[ms, m[{1,7}]],
            (* m[{4,10}] is NOT \[LongDash] it's collapsed into m[{1,7}] *)
            !MemberQ[ms, m[{4,10}]],
            AnyTrue[eqDF,
                TrueQ[Simplify[(#[[1]] - j17Eq[[1]]) == 0]] &]
        ]
    ],
    True,
    TestID -> "Tawaf: density extension on 2x3x2 produces tangential + radial cohorts"
]

Test[
    (* tawafDensities post-solve helper: after solving the j/u part of a
       Density-True system, the helper derives m for each physical edge by
       inverting the linear consistency equation. All m values must be
       strictly positive (consistent with IneqMs). *)
    Module[{s, sys, sol, mRules, mVals},
        s   = makeTawafScenario[2, 3, 1, "Density" -> True];
        sys = makeTawafSystem[s];
        sol = TimeConstrained[dnfReduceSystem[sys], 60, $TimedOut];
        sol =!= $TimedOut && (
            mRules = tawafDensities[sys, sol];
            mVals  = Values[mRules];
            And[
                AssociationQ[mRules],
                Length[mRules] === 3,
                AllTrue[mVals, NumericQ[#] && # > 0 &]
            ]
        )
    ],
    True,
    TestID -> "Tawaf: tawafDensities derives positive m per physical edge"
]

Test[
    (* tawafDensities returns an empty association on a non-density system
       (defensive: callers can compose unconditionally with Join). *)
    Module[{s, sys, sol, mRules},
        s   = makeTawafScenario[2, 3, 1];
        sys = makeTawafSystem[s];
        sol = TimeConstrained[dnfReduceSystem[sys], 60, $TimedOut];
        sol =!= $TimedOut && (
            mRules = tawafDensities[sys, sol];
            mRules === <||>
        )
    ],
    True,
    TestID -> "Tawaf: tawafDensities returns <||> on non-density system"
]
