(* Wolfram Language Test file *)
(* Tests for package loading and context hygiene *)

(* Test: Key public symbols exist in MFGraphs` context *)
Test[
    NameQ["MFGraphs`DataToEquations"] && NameQ["MFGraphs`CriticalCongestionSolver"] &&
    NameQ["MFGraphs`SolveMFG"] &&
    NameQ["MFGraphs`GetExampleData"] && NameQ["MFGraphs`DNFReduce"] &&
    NameQ["MFGraphs`GetKirchhoffLinearSystem"] &&
    NameQ["MFGraphs`makeScenario"] && NameQ["MFGraphs`scenarioQ"]
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
    DataG === GetExampleData && Data2Equations === DataToEquations && FinalStep === DNFSolveStep &&
    RemoveDuplicates === DeduplicateByComplexity && ReplaceSolution === SubstituteSolution
    ,
    True
    ,
    TestID -> "Package loading: backward compatibility aliases"
]
