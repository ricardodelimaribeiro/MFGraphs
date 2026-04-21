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
    Module[{data, d2e, result, unresolved},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        d2e = Join[d2e, <|"CriticalNumericBackendMode" -> False|>];
        result = Quiet[CriticalCongestionSolver[d2e, "SymbolicTimeLimit" -> 120.]];
        unresolved = Lookup[result, "UnresolvedEquations", Missing["NotAvailable"]];
        AssociationQ[result] &&
        TrueQ[IsFeasible[result]] &&
        AssociationQ[Lookup[result, "AssoCritical", Missing["NotAvailable"]]] &&
        unresolved =!= True &&
        unresolved =!= None &&
        unresolved =!= Missing["NotAvailable"]
    ]
    ,
    True
    ,
    TestID -> "Critical symbolic path: underdetermined case exports unresolved equations"
]
