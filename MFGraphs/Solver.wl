(* Wolfram Language package *)
(*
   Solver.wl: symbolic solver for mfgSystem objects.

   Current scope: ReduceSystem collects the structural equations,
   flow-balance constraints, non-negativity inequalities, and
   complementarity alternatives from an mfgSystem, substitutes the
   exit-value boundary conditions, then calls Reduce over the Reals.
   Switching-cost inequalities are not included.
*)

BeginPackage["MFGraphs`"];

ReduceSystem::usage =
"ReduceSystem[sys] reduces the structural equations, flow balance, \
non-negativity constraints, and complementarity conditions of the \
mfgSystem sys using Reduce over the Reals. Switching-cost inequalities \
(IneqSwitchingByVertex, AltOptCond) are not included in this pass.";

Begin["`Private`"];

ReduceSystem[sys_?mfgSystemQ] :=
    Module[{eqEntryIn, eqSplit, eqGather, eqHJ,
            ineqJs, ineqJts, altFlows, altTrans,
            ruleExitVals, constraints, allVars},

        eqEntryIn    = SystemData[sys, "EqEntryIn"];
        eqSplit      = SystemData[sys, "EqBalanceSplittingFlows"];
        eqGather     = SystemData[sys, "EqBalanceGatheringFlows"];
        eqHJ         = SystemData[sys, "EqGeneral"];
        ineqJs       = SystemData[sys, "IneqJs"];
        ineqJts      = SystemData[sys, "IneqJts"];
        altFlows     = SystemData[sys, "AltFlows"];
        altTrans     = SystemData[sys, "AltTransitionFlows"];
        ruleExitVals = Normal @ SystemData[sys, "RuleExitValues"];

        constraints = And[
            And @@ eqEntryIn,
            eqSplit, eqGather, eqHJ,
            ineqJs, ineqJts,
            altFlows, altTrans
        ] /. ruleExitVals;

        allVars = Select[
            Variables[constraints /. {Equal -> List, Or -> List, And -> List}],
            MatchQ[#, _j | _u] &
        ];

        Reduce[constraints, allVars, Reals]
    ];

End[];

EndPackage[];
