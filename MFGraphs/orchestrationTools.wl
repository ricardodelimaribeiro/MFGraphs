(* Wolfram Language package *)
(* orchestrationTools.wl: High-level wrappers connecting scenario kernels to solvers. *)

BeginPackage["orchestrationTools`", {
    "primitives`", 
    "scenarioTools`", 
    "unknownsTools`", 
    "systemTools`", 
    "solversTools`"
}];

solveScenario::usage =
"solveScenario[s] automatically constructs exact symbolic unknowns, builds the \
structural system, and calls dnfReduceSystem. \
solveScenario[{s1, s2, ...}] solves multiple populations (scenarios) and \
returns a list of solutions. \
solveScenario[..., solver] uses the specified solver function (e.g., reduceSystem).";

SolveMFG::usage =
"SolveMFG[s] solves a typed scenario object by delegating to solveScenario. \
SolveMFG[assoc] provides backward compatibility for legacy raw-association \
solving. It constructs a scenario and delegates to solveScenario. \
SolveMFG[s, solver] or SolveMFG[assoc, solver] uses the specified solver function.";

Begin["`Private`"];

(* Single scenario orchestration *)
solveScenario[s_?scenarioQ, solver_:dnfReduceSystem] :=
    Module[{unk, sys},
        unk = makeSymbolicUnknowns[s];
        sys = makeSystem[s, unk];
        solver[sys]
    ];

(* Multi-population orchestration *)
solveScenario[scenarios_List, solver_:dnfReduceSystem] :=
    solveScenario[#, solver] & /@ scenarios;

SolveMFG[s_?scenarioQ, solver_:dnfReduceSystem] :=
    solveScenario[s, solver];

SolveMFG[assoc_Association, solver_:dnfReduceSystem] :=
    Module[{s},
        s = makeScenario[assoc];
        If[FailureQ[s],
            Message[SolveMFG::invalid, s["Message"]];
            Return[s]
        ];
        solveScenario[s, solver]
    ];

SolveMFG::invalid = "Legacy SolveMFG input failed scenario conversion: `1`";

End[];

EndPackage[];
