(* Wolfram Language package *)
(* Unknowns: construction helpers for MFGraphs symbolic unknown families. *)

BeginPackage["MFGraphs`"];

unknowns::usage = "unknowns[assoc] is the typed head for MFGraphs unknown-variable bundles.";

unknownsQ::usage = "unknownsQ[x] returns True iff x is an unknowns[assoc_Association] object.";

UnknownsData::usage = "UnknownsData[u, key] returns key from unknowns object u; UnknownsData[u] returns the underlying association.";

makeUnknowns::usage = "makeUnknowns[s] returns unknowns[<|\"js\" -> ..., \"jts\" -> ..., \"us\" -> ...|>] for scenario s. \"js\" are flow unknowns j[v,e], \"jts\" are transition-flow unknowns j[v,e1,e2], and \"us\" are value-function unknowns u[v,e].";

Begin["`Private`"];

iCanonicalUPair[pair : {a_, b_}] :=
    Module[{bName = Quiet @ Check[SymbolName[b], ""]},
        If[StringStartsQ[bName, "auxEntry"] || StringStartsQ[bName, "auxExit"],
            {b, a},
            pair
        ]
    ];

MakeUnknownsFromPairsTriples[auxPairs_List, auxTriples_List] :=
    unknowns[
        <|
            "js" -> (j[Sequence @@ #] & /@ auxPairs),
            "jts" -> (j[Sequence @@ #] & /@ auxTriples),
            "us" -> (u[Sequence @@ #] & /@ (iCanonicalUPair /@ auxPairs)),
            "auxPairs" -> auxPairs,
            "auxTriples" -> auxTriples
        |>
    ];

unknownsQ[unknowns[_Association]] := True;
unknownsQ[_] := False;

UnknownsData[unknowns[assoc_Association]] := assoc;
UnknownsData[unknowns[assoc_Association], key_] := Lookup[assoc, key, Missing["KeyAbsent", key]];

makeUnknowns[s_] :=
    Module[{model, topology, rawAssoc, auxPairs, auxTriples},
        rawAssoc = Which[
            MFGraphs`scenarioQ[s], ScenarioData[s],
            MatchQ[s, _[ _Association]], First[s],
            AssociationQ[s], s,
            True, $Failed
        ];

        If[!AssociationQ[rawAssoc],
            Return[Failure["makeUnknowns", <|"Message" -> "Input must be a scenario or a model association."|>], Module]
        ];

        topology = If[MFGraphs`scenarioQ[s], 
            ScenarioData[s, "Topology"], 
            Lookup[rawAssoc, "Topology", $Failed]
        ];
        
        If[!AssociationQ[topology],
            model = Lookup[rawAssoc, "Model", rawAssoc];
            If[!AssociationQ[model],
                Return[Failure["makeUnknowns", <|"Message" -> "Input must be a scenario or a model association."|>], Module]
            ];
            topology = MFGraphs`BuildAuxiliaryTopology[model];
        ];

        If[topology === $Failed,
            Return[Failure["makeUnknowns", <|"Message" -> "Could not build auxiliary topology."|>], Module]
        ];

        (* Extract pairs and triples from pre-calculated topology *)
        auxPairs = Lookup[topology, "AuxPairs", MFGraphs`DeriveAuxPairs[topology]];
        auxTriples = Lookup[topology, "AuxTriples", MFGraphs`BuildAuxTriples[topology["AuxiliaryGraph"]]];
        
        MakeUnknownsFromPairsTriples[auxPairs, auxTriples]
    ];

End[];

EndPackage[];
