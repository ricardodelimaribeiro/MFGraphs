(* Wolfram Language Test file *)
(* Tests for package loading and context hygiene *)

(* Test: Key public symbols exist in MFGraphs` context *)
Test[
    NameQ["MFGraphs`DataToEquations"] && NameQ["MFGraphs`CriticalCongestionSolver"] &&
    NameQ["MFGraphs`NonLinearSolver"] && NameQ["MFGraphs`MonotoneSolverFromData"] &&
    NameQ["MFGraphs`SolveMFG"] &&
    NameQ["MFGraphs`GetExampleData"] && NameQ["MFGraphs`DNFReduce"] &&
    NameQ["MFGraphs`GetKirchhoffLinearSystem"] && NameQ["MFGraphs`V"] &&
    NameQ["MFGraphs`alpha"] && NameQ["MFGraphs`g"] &&
    NameQ["MFGraphs`WithHamiltonianFunctions"] &&
    NameQ["MFGraphs`makeScenario"] && NameQ["MFGraphs`scenarioQ"]
    ,
    True
    ,
    TestID -> "Package loading: public symbols exist"
]

(* Test: Private symbols are NOT in MFGraphs` context *)
Test[
    !NameQ["MFGraphs`$SolveCache"] && !NameQ["MFGraphs`$ReduceCache"] &&
    !NameQ["MFGraphs`sortByComplexity"] && !NameQ["MFGraphs`TransitionsAt"] &&
    !NameQ["MFGraphs`BuildMonotoneField"] && !NameQ["MFGraphs`MakeSolverResult"]
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

Test[
    {
        V[0, UndirectedEdge[1, 2]],
        alpha[UndirectedEdge[1, 2]],
        g[2, UndirectedEdge[1, 2]]
    } ===
    {0, 1, -1/4}
    ,
    True
    ,
    TestID -> "Package loading: default Hamiltonian hooks"
]

Test[
    WithHamiltonianFunctions[
        Function[{x, edge}, x + 3],
        Function[edge, 2],
        Function[{m, edge}, -3 m],
        {
            V[2, UndirectedEdge[1, 2]],
            alpha[UndirectedEdge[1, 2]],
            g[2, UndirectedEdge[1, 2]]
        }
    ] === {5, 2, -6}
    ,
    True
    ,
    TestID -> "Package loading: temporary Hamiltonian overrides"
]
