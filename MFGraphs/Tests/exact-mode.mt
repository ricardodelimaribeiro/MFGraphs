(* Wolfram Language Test file *)

If[!MemberQ[$Packages, "MFGraphs`"],
    Get["/Users/ribeirrd/Documents/GitHub/MFGraphs/MFGraphs/MFGraphs.wl"]
];

Test[
    Module[{data, d2e, result, asso, region},
        data = GetExampleData[12] /. {I1 -> 100, U1 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[CriticalCongestionSolver[d2e, "ExactMode" -> True]];
        asso = Lookup[result, "AssoCritical", Missing["NotAvailable"]];
        region = Lookup[result, "SymbolicRegion", Missing["NotAvailable"]];
        AssociationQ[result] &&
        TrueQ[IsFeasible[result]] &&
        AssociationQ[asso] &&
        Length[asso] > 0 &&
        region =!= None &&
        !MissingQ[region] &&
        Lookup[result, "NumericBackendUsed", Missing["NotAvailable"]] === False &&
        Lookup[result, "JFirstBackendUsed", Missing["NotAvailable"]] === False
    ]
    ,
    True
    ,
    TestID -> "ExactMode: case 12 returns symbolic feasible envelope without numeric backends"
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
        Lookup[result, "NumericBackendUsed", Missing["NotAvailable"]] === False &&
        Lookup[result, "JFirstBackendUsed", Missing["NotAvailable"]] === False
    ]
    ,
    True
    ,
    TestID -> "ExactMode: underdetermined case exposes SymbolicRegion"
]
