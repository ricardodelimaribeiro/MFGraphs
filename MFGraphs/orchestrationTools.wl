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
"solveScenario[s] automatically constructs unknowns, builds the structural \
system, and calls reduceSystem. \
solveScenario[{s1, s2, ...}] solves multiple populations (scenarios) and \
returns a list of solutions. \
solveScenario[..., solver] uses the specified solver function (e.g., dnfReduceSystem).";

SolveMFG::usage =
"SolveMFG[assoc] provides backward compatibility for legacy raw-association \
solving. It constructs a scenario and delegates to solveScenario.";

Begin["`Private`"];

(* Single scenario orchestration *)
solveScenario[s_?scenarioQ, solver_:reduceSystem] :=
    Module[{unk, sys},
        unk = makeUnknowns[s];
        sys = makeSystem[s, unk];
        solver[sys]
    ];

(* Multi-population orchestration *)
solveScenario[scenarios_List, solver_:reduceSystem] :=
    solveScenario[#, solver] & /@ scenarios;

SolveMFG[assoc_Association] :=
    Module[{s},
        s = makeScenario[assoc];
        If[FailureQ[s],
            Message[SolveMFG::invalid, s["Message"]];
            Return[s]
        ];
        solveScenario[s]
    ];

SolveMFG::invalid = "Legacy SolveMFG input failed scenario conversion: `1`";

End[];

EndPackage[];
