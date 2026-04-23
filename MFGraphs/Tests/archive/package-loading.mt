(* Wolfram Language Test file *)
(* Tests for package loading and context hygiene *)

(* Test: Key public symbols exist in MFGraphs` context *)
Test[
    NameQ["MFGraphs`DataToEquations"] && NameQ["MFGraphs`CriticalCongestionSolver"] &&
    NameQ["MFGraphs`SolveMFG"] &&
    NameQ["MFGraphs`GetExampleData"] && NameQ["MFGraphs`DNFReduce"] &&
    NameQ["MFGraphs`GetKirchhoffLinearSystem"] &&
    NameQ["MFGraphs`makeScenario"] && NameQ["MFGraphs`scenarioQ"] &&
    NameQ["MFGraphs`ScenarioByKey"]
    ,
    True
    ,
    TestID -> "Package loading: public symbols exist"
]

(* Test: Private symbols are NOT in MFGraphs` context *)
Test[
    !NameQ["MFGraphs`$SolveCache"] && !NameQ["MFGraphs`$ReduceCache"] &&
    !NameQ["MFGraphs`sortByComplexity"] && !NameQ["MFGraphs`TransitionsAt"]
    ,
    True
    ,
    TestID -> "Package loading: private symbols not leaked"
]

(* Test: Backward compatibility aliases resolve correctly *)
Test[
    DataG === GetExampleData && FinalStep === DNFSolveStep &&
    RemoveDuplicates === DeduplicateByComplexity && ReplaceSolution === SubstituteSolution
    ,
    True
    ,
    TestID -> "Package loading: backward compatibility aliases"
]

(* Test: Deprecated Data2Equations remains behavior-compatible with DataToEquations *)
Test[
    Module[{data, d2e1, d2e2},
        data = GetExampleData[3] /. {I1 -> 2, U1 -> 0};
        d2e1 = DataToEquations[data];
        d2e2 = Quiet[Data2Equations[data], Data2Equations::deprecated];
        AssociationQ[d2e1] && d2e1 === d2e2
    ]
    ,
    True
    ,
    TestID -> "Package loading: Data2Equations deprecated wrapper is behavior-compatible"
]

(* Test: Data2Equations emits deprecation message *)
Test[
    Module[{data},
        data = GetExampleData[3] /. {I1 -> 2, U1 -> 0};
        Quiet[
            Check[
                Data2Equations[data];
                False,
                True,
                Data2Equations::deprecated
            ],
            Data2Equations::deprecated
        ]
    ]
    ,
    True
    ,
    TestID -> "Package loading: Data2Equations emits deprecation warning"
]
