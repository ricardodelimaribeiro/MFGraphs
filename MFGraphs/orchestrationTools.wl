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
solveScenario[..., solver] uses the specified solver function (e.g., reduceSystem). \
Results are memoized per (scenario, solver) within the session; call \
clearSolveCache[] to wipe.";

SolveMFG::usage =
"SolveMFG[s] solves a typed scenario object by delegating to solveScenario. \
SolveMFG[assoc] provides backward compatibility for legacy raw-association \
solving. It constructs a scenario and delegates to solveScenario. \
SolveMFG[s, solver] or SolveMFG[assoc, solver] uses the specified solver function.";

clearSolveCache::usage =
"clearSolveCache[] empties solveScenario's session-scoped memoization cache. \
Call between benchmark passes to measure cold-start cost, or to release memory \
in long-running sessions.";

Begin["`Private`"];

(* Session-scoped memoization. Solvers in this package (dnfReduceSystem,
   optimizedDNFReduceSystem, activeSetReduceSystem, ...) are pure functions of
   their mfgSystem input, so caching solveScenario[s, solver] by Hash of the
   scenario payload + solver symbol is sound. The first call pays the full
   makeSystem + solve cost; subsequent calls on the same (s, solver) return in
   O(hash). Direct calls into dnfReduceSystem[sys] from benchmarks bypass this
   cache, as intended. *)
$solveScenarioCache = <||>;

clearSolveCache[] := ($solveScenarioCache = <||>;);

(* Single scenario orchestration *)
solveScenario[s_?scenarioQ, solver_:dnfReduceSystem] :=
    Module[{key, cached, unk, sys, result},
        key = Hash[{scenarioData[s], solver}];
        cached = Lookup[$solveScenarioCache, key];
        If[!MissingQ[cached], Return[cached, Module]];
        unk = makeSymbolicUnknowns[s];
        sys = makeSystem[s, unk];
        result = solver[sys];
        $solveScenarioCache[key] = result;
        result
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
