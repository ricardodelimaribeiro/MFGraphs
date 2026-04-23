(* Wolfram Language package *)
(* Unknowns: construction helpers for MFGraphs symbolic unknown families. *)

unknowns::usage = "unknowns[assoc] is the typed head for MFGraphs unknown-variable bundles.";

unknownsQ::usage = "unknownsQ[x] returns True iff x is an unknowns[assoc_Association] object.";

UnknownsData::usage = "UnknownsData[u, key] returns key from unknowns object u; UnknownsData[u] returns the underlying association.";

makeUnknowns::usage = "makeUnknowns[s] returns unknowns[<|\"js\" -> ..., \"jts\" -> ..., \"us\" -> ...|>] for scenario s. \"js\" are flow unknowns j[v,e], \"jts\" are transition-flow unknowns j[v,e1,e2], and \"us\" are value-function unknowns u[v,e].";

Begin["`Private`"];

DeriveAuxPairs[topology_Association] :=
    Module[{graph, halfPairs, inAuxEntryPairs, outAuxExitPairs, inAuxExitPairs, 
            outAuxEntryPairs, pairs},
        graph = topology["Graph"];
        halfPairs = List @@@ EdgeList[graph];
        inAuxEntryPairs = List @@@ topology["AuxEntryEdges"];
        outAuxExitPairs = List @@@ topology["AuxExitEdges"];
        inAuxExitPairs = Reverse /@ outAuxExitPairs;
        outAuxEntryPairs = Reverse /@ inAuxEntryPairs;
        pairs = Join[halfPairs, Reverse /@ halfPairs];
        Join[inAuxEntryPairs, outAuxEntryPairs, inAuxExitPairs, outAuxExitPairs, pairs]
    ];

MakeUnknownsFromPairsTriples[auxPairs_List, auxTriples_List] :=
    unknowns[
        <|
            "js" -> (j[Sequence @@ #] & /@ auxPairs),
            "jts" -> (j[Sequence @@ #] & /@ auxTriples),
            "us" -> (u[Sequence @@ #] & /@ auxPairs),
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

        (* Combinatorial creation of pairs and triples moved here *)
        auxPairs = DeriveAuxPairs[topology];
        auxTriples = If[
            MFGraphs`scenarioQ[s] && AssociationQ[topology] && KeyExistsQ[topology, "AuxTriples"],
            topology["AuxTriples"],
            Lookup[topology, "AuxTriples", MFGraphs`BuildAuxTriples[topology["AuxiliaryGraph"]]]
        ];
        
        MakeUnknownsFromPairsTriples[auxPairs, auxTriples]
    ];

End[];
