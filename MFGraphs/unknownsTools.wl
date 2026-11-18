(* Wolfram Language package *)
(* unknownsTools.wl — construction helpers for MFGraphs symbolic unknown families. *)

BeginPackage["unknownsTools`", {"primitives`", "utilities`", "scenarioTools`"}]

symbolicUnknowns::usage =
"symbolicUnknowns[assoc] is the typed head for topology-derived exact symbolic \
variable bundles used by MFGraphs structural graph systems.";

symbolicUnknownsQ::usage =
"symbolicUnknownsQ[x] returns True iff x is a symbolicUnknowns[assoc_Association] object.";

symbolicUnknownsData::usage =
"symbolicUnknownsData[u, key] returns key from symbolicUnknowns object u; \
symbolicUnknownsData[u] returns the underlying association.";

makeSymbolicUnknowns::usage =
"makeSymbolicUnknowns[s] returns symbolicUnknowns[<|\"Js\" -> ..., \"Jts\" -> ..., \
\"Us\" -> ...|>] for scenario s. \"Js\" are exact flow variables j[a,b], \
\"Jts\" are exact transition-flow variables j[r,i,w], and \"Us\" are exact \
value-function variables u[a,b]. A scenario with cached topology is used directly. \
Any other input (a scenario without topology, a raw scenario association, or a \
bare model association) is first passed through makeScenario, so it is validated \
the same way, and makeScenario's Failure is returned when it is invalid.";

unknown::usage =
"unknown is reserved for future numeric MFGraphs solver fields over grids or \
layouts; it is not used by the exact symbolic graph-system constructors.";

unknowns::usage =
"unknowns is reserved for future numeric collections of unknown fields over \
grids or layouts, following the Maydan-style numeric solver boundary. Use \
symbolicUnknowns for current exact symbolic graph systems.";

Begin["`Private`"]

canonicalUPair::usage =
"canonicalUPair[{a,b}] returns the canonical pair used for value-function \
variables, orienting auxiliary entry/exit pairs consistently.";

makeSymbolicUnknownsFromPairsTriples::usage =
"makeSymbolicUnknownsFromPairsTriples[auxPairs, auxTriples] builds the typed \
symbolicUnknowns object from auxiliary pairs and triples.";

canonicalUPair[pair : {a_, b_}] :=
    If[MatchQ[b, _auxEntry | _auxExit],
        {b, a},
        pair
    ];

makeSymbolicUnknownsFromPairsTriples[auxPairs_List, auxTriples_List] :=
    symbolicUnknowns[
        <|
            "Js"        -> (j[Sequence @@ #] & /@ auxPairs),
            "Jts"       -> (j[Sequence @@ #] & /@ auxTriples),
            "Us"        -> (u[Sequence @@ #] & /@ (canonicalUPair /@ auxPairs)),
            "AuxPairs"  -> auxPairs,
            "AuxTriples" -> auxTriples
        |>
    ];

symbolicUnknownsQ[x_] := mfgTypedQ[x, symbolicUnknowns];

symbolicUnknownsData[u_] := mfgData[u];
symbolicUnknownsData[u_, key_] := mfgData[u, key];

makeSymbolicUnknowns[s_] :=
    Module[{topology, rawAssoc, built},
        topology = If[scenarioQ[s], scenarioData[s, "Topology"], Missing["NotAScenario"]];

        (* Anything other than a scenario with cached topology goes through
           makeScenario, so this path cannot accept a model that makeScenario
           would reject. *)
        If[!AssociationQ[topology],
            rawAssoc = Which[
                scenarioQ[s], scenarioData[s],
                MatchQ[s, _[_Association]], First[s],
                AssociationQ[s], s,
                True, $Failed
            ];
            If[!AssociationQ[rawAssoc],
                Return[Failure["makeSymbolicUnknowns", <|"Message" -> "Input must be a scenario or a model association."|>], Module]
            ];
            built = makeScenario[If[KeyExistsQ[rawAssoc, "Model"], rawAssoc, <|"Model" -> rawAssoc|>]];
            If[FailureQ[built], Return[built, Module]];
            topology = scenarioData[built, "Topology"]
        ];

        makeSymbolicUnknownsFromPairsTriples[topology["AuxPairs"], topology["AuxTriples"]]
    ];

End[]

EndPackage[]
