(* Wolfram Language Test file *)

If[!MemberQ[$Packages, "MFGraphs`"],
    Get[
        FileNameJoin[
            {
                DirectoryName[$InputFileName],
                "..",
                "MFGraphs.wl"
            }
        ]
    ]
];

Test[
    Module[{data, d2e, result, region, asso},
        data = GetExampleData[12] /. {I1 -> 100, U1 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[CriticalCongestionSolver[d2e, "ExactMode" -> True]];
        region = Lookup[result, "SymbolicRegion", Missing["NotAvailable"]];
        asso = Lookup[result, "AssoCritical", Missing["NotAvailable"]];
        AssociationQ[result] &&
        TrueQ[IsFeasible[result]] &&
        AssociationQ[asso] &&
        Length[asso] > 0 &&
        !MissingQ[region] &&
        Lookup[result, "NumericBackendUsed", Missing["NotAvailable"]] === False &&
        Lookup[result, "JFirstBackendUsed", Missing["NotAvailable"]] === False
    ]
    ,
    True
    ,
    TestID -> "ExactMode: fully determined case yields no SymbolicRegion"
]

Test[
    Module[{data, d2e, result, region},
        data = GetExampleData["triangle with two exits"] /. {I1 -> 100, U1 -> 0, S1 -> 1};
        d2e = DataToEquations[data];
        result = Quiet[CriticalCongestionSolver[d2e, "ExactMode" -> True]];
        region = Lookup[result, "SymbolicRegion", Missing["NotAvailable"]];
        AssociationQ[result] &&
        TrueQ[IsFeasible[result]] &&
        region =!= None &&
        !MissingQ[region] &&
        Lookup[result, "NumericBackendUsed", Missing["NotAvailable"]] === False
    ]
    ,
    True
    ,
    TestID -> "ExactMode: underdetermined case exposes SymbolicRegion"
]
