(* ::Package:: *)

(*
  This script compares the solutions of CriticalCongestionSolver, NonLinearSolver,
  and MonotoneSolver for a critical congestion problem.

  All three solvers should find the same solution when their parameters are set to match
  the assumptions of the critical congestion model (V=0, alpha=1, g=-1/m^2).
*)

(* --- Setup: Load MFGraphs --- *)
PrependTo[$Path, Directory[]];
Needs["MFGraphs`"];
MFGraphs`$MFGraphsVerbose = True;


(* --- 1. Define the problem --- *)
Print["1. Setting up the problem from example data..."];
data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
d2e = DataToEquations[data];
Print["Problem setup complete."];
Print[""];

(* --- 2. Solve with CriticalCongestionSolver --- *)
Print["2. Solving with CriticalCongestionSolver..."];
criticalResult = CriticalCongestionSolver[d2e];
criticalSolution = criticalResult["AssoCritical"];
If[!AssociationQ[criticalSolution],
    Print["CriticalCongestionSolver failed to find a solution."];
    Exit[1];
];
Print["CriticalCongestionSolver found a solution."];
Print[""];


(* --- 3. Solve with NonLinearSolver --- *)
Print["3. Solving with NonLinearSolver..."];
(*
  For a critical congestion problem, we use the default Hamiltonian:
  - Potential V(x, e) = 0
  - Congestion exponent alpha(e) = 1
  - Interaction g(m, e) = -1/m^2
*)
potentialFunction = Function[{x, edge}, 0];
congestionExponentFunction = Function[edge, 1];
interactionFunction = Function[{m, edge}, -1/m^2];

nonlinearResult = NonLinearSolver[d2e,
    "PotentialFunction" -> potentialFunction,
    "CongestionExponentFunction" -> congestionExponentFunction,
    "InteractionFunction" -> interactionFunction
];
nonlinearSolution = nonlinearResult["AssoNonCritical"];

If[!AssociationQ[nonlinearSolution],
    Print["NonLinearSolver failed to find a solution."];
    Exit[1];
];
Print["NonLinearSolver found a solution."];
Print[""];


(* --- 4. Solve with MonotoneSolver --- *)
Print["4. Solving with MonotoneSolver..."];
monotoneResult = MonotoneSolver[d2e,
    "PotentialFunction" -> potentialFunction,
    "CongestionExponentFunction" -> congestionExponentFunction,
    "InteractionFunction" -> interactionFunction,
    "TimeSteps" -> 100
];

If[!AssociationQ[monotoneResult],
    Print["MonotoneSolver failed to find a solution."];
    Exit[1];
];
Print["MonotoneSolver found a solution."];
Print[""];


(* --- 5. Compare the solutions --- *)
Print["5. Comparing solutions..."];
tolerance = 10^-6;
allPass = True;

(* Helper: compare two solution associations *)
comparePair[name1_, sol1_, name2_, sol2_] :=
    Module[{keys, vals1, vals2, diff},
        keys = Sort[Intersection[Keys[sol1], Keys[sol2]]];
        vals1 = keys /. sol1;
        vals2 = keys /. sol2;
        diff = Norm[vals1 - vals2];
        Print["  ", name1, " vs ", name2, ": norm = ", diff];
        If[diff >= tolerance,
            Print["  FAIL: difference exceeds tolerance (", tolerance, ")"];
            allPass = False;
        ];
    ];

comparePair["Critical", criticalSolution, "NonLinear", nonlinearSolution];
comparePair["Critical", criticalSolution, "Monotone", monotoneResult];
comparePair["NonLinear", nonlinearSolution, "Monotone", monotoneResult];
Print[""];

If[allPass,
    Print["SUCCESS: All solver pairs are consistent within tolerance (", tolerance, ")."];
    Exit[0];
,
    Print["FAILURE: Some solver pairs are NOT consistent."];
    Exit[1];
]
