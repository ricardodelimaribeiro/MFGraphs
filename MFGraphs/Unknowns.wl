(* Wolfram Language package *)
(* Unknowns: construction helpers for MFGraphs symbolic unknown families. *)

unknowns::usage = "unknowns[assoc] is the typed head for MFGraphs unknown-variable bundles.";

unknownsQ::usage = "unknownsQ[x] returns True iff x is an unknowns[assoc_Association] object.";

UnknownsData::usage = "UnknownsData[u, key] returns key from unknowns object u; UnknownsData[u] returns the underlying association.";

makeUnknowns::usage = "makeUnknowns[s] returns unknowns[<|\"js\" -> ..., \"jts\" -> ..., \"us\" -> ...|>] for scenario s. \"js\" are flow unknowns j[v,e], \"jts\" are transition-flow unknowns j[v,e1,e2], and \"us\" are value-function unknowns u[v,e].";

Begin["`Private`"];

MakeUnknownsFromPairsTriples[auxPairs_List, auxTriples_List] :=
    unknowns[
        <|
            "js" -> (j[Sequence @@ #] & /@ auxPairs),
            "jts" -> (j[Sequence @@ #] & /@ auxTriples),
            "us" -> (u[Sequence @@ #] & /@ auxPairs)
        |>
    ];

unknownsQ[unknowns[_Association]] := True;
unknownsQ[_] := False;

UnknownsData[unknowns[assoc_Association]] := assoc;
UnknownsData[unknowns[assoc_Association], key_] := Lookup[assoc, key, Missing["KeyAbsent", key]];

makeUnknowns[s_] :=
    Module[{model, topology},
        model = If[scenarioQ[s],
            ScenarioData[s, "Model"],
            If[AssociationQ[s], Lookup[s, "Model", s], $Failed]
        ];
        
        If[!AssociationQ[model],
            Return[
                Failure["makeUnknowns",
                    <|"Message" -> "Input must be a scenario or a model association."|>
                ],
                Module
            ]
        ];
        
        topology = BuildAuxiliaryTopology[model];
        If[topology === $Failed,
            Return[
                Failure["makeUnknowns",
                    <|"Message" -> "Could not build auxiliary topology from model."|>
                ],
                Module
            ]
        ];
        
        MakeUnknownsFromPairsTriples[topology["AuxPairs"], topology["AuxTriples"]]
    ];

End[];
