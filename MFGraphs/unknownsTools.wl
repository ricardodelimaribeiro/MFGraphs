(* Wolfram Language package *)
(* unknownsTools.wl — construction helpers for MFGraphs symbolic unknown families. *)

BeginPackage["unknownsTools`", {"primitives`", "scenarioTools`"}]

unknowns::usage = "unknowns[assoc] is the typed head for MFGraphs unknown-variable bundles.";

unknownsQ::usage = "unknownsQ[x] returns True iff x is an unknowns[assoc_Association] object.";

unknownsData::usage = "unknownsData[u, key] returns key from unknowns object u; unknownsData[u] returns the underlying association.";

makeUnknowns::usage = "makeUnknowns[s] returns unknowns[<|\"Js\" -> ..., \"Jts\" -> ..., \"Us\" -> ...|>] for scenario s. \"Js\" are flow unknowns j[v,e], \"Jts\" are transition-flow unknowns j[v,e1,e2], and \"Us\" are value-function unknowns u[v,e].";

Begin["`Private`"]

canonicalUPair[pair : {a_, b_}] :=
    If[StringQ[b] && (StringStartsQ[b, "auxEntry"] || StringStartsQ[b, "auxExit"]),
        {b, a},
        pair
    ];

makeUnknownsFromPairsTriples[auxPairs_List, auxTriples_List] :=
    unknowns[
        <|
            "Js"        -> (j[Sequence @@ #] & /@ auxPairs),
            "Jts"       -> (j[Sequence @@ #] & /@ auxTriples),
            "Us"        -> (u[Sequence @@ #] & /@ (canonicalUPair /@ auxPairs)),
            "AuxPairs"  -> auxPairs,
            "AuxTriples" -> auxTriples
        |>
    ];

unknownsQ[unknowns[_Association]] := True;
unknownsQ[_] := False;

unknownsData[unknowns[assoc_Association]] := assoc;
unknownsData[unknowns[assoc_Association], key_] := Lookup[assoc, key, Missing["KeyAbsent", key]];

makeUnknowns[s_] :=
    Module[{model, topology, rawAssoc, auxPairs, auxTriples},
        rawAssoc = Which[
            scenarioQ[s], scenarioData[s],
            MatchQ[s, _[ _Association]], First[s],
            AssociationQ[s], s,
            True, $Failed
        ];

        If[!AssociationQ[rawAssoc],
            Return[Failure["makeUnknowns", <|"Message" -> "Input must be a scenario or a model association."|>], Module]
        ];

        topology = If[scenarioQ[s],
            scenarioData[s, "Topology"],
            Lookup[rawAssoc, "Topology", $Failed]
        ];

        If[!AssociationQ[topology],
            model = Lookup[rawAssoc, "Model", rawAssoc];
            If[!AssociationQ[model],
                Return[Failure["makeUnknowns", <|"Message" -> "Input must be a scenario or a model association."|>], Module]
            ];
            topology = buildAuxiliaryTopology[model];
        ];

        If[topology === $Failed,
            Return[Failure["makeUnknowns", <|"Message" -> "Could not build auxiliary topology."|>], Module]
        ];

        auxPairs   = Lookup[topology, "AuxPairs",   deriveAuxPairs[topology]];
        auxTriples = Lookup[topology, "AuxTriples",  buildAuxTriples[topology["AuxiliaryGraph"]]];

        makeUnknownsFromPairsTriples[auxPairs, auxTriples]
    ];

End[]

EndPackage[]
