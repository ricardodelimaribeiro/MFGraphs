(* Tests for the opt-in numericOracle subpackage and its high-level wrapper
   solveScenarioWithOracle. numericOracle is NOT loaded by Needs["MFGraphs`"],
   so this suite begins with an explicit Needs["numericOracle`"]. *)

Needs["MFGraphs`"];
Needs["numericOracle`"];

Test[
    NameQ["numericOracle`numericOracleClassify"] &&
    NameQ["numericOracle`solveScenarioWithOracle"],
    True,
    TestID -> "numericOracle: public symbols exist"
]

(* Happy path: linear 3-chain. Single entry at vertex 1, single exit at vertex 3.
   Fully determined; the active-set solver returns a list of rules whether or
   not the oracle prunes anything, and the wrapper must return that list. *)
Test[
    Module[{s, result},
        s = getExampleScenario[3, {{1, 1}}, {{3, 0}}];
        result = solveScenarioWithOracle[s];
        ListQ[result] && AllTrue[result, MatchQ[#, _Rule | _RuleDelayed] &]
    ],
    True,
    TestID -> "solveScenarioWithOracle: returns rule list for fully-determined 3-chain"
]

(* Wrapper-vs-direct equivalence: the oracle path must produce the same
   solution (modulo rule order) as the direct activeSetReduceSystem path on
   a scenario where the oracle does not over-prune. *)
Test[
    Module[{s, unk, sys, sym, direct, oracleResult,
            canon},
        s    = getExampleScenario[3, {{1, 1}}, {{3, 0}}];
        unk  = makeSymbolicUnknowns[s];
        sys  = makeSystem[s, unk];
        sym  = addSymmetryEqualities[sys, s];
        direct       = activeSetReduceSystem[sym];
        oracleResult = solveScenarioWithOracle[s];
        canon[r_] := If[ListQ[r], SortBy[r, ToString[First[#]] &], r];
        canon[direct] === canon[oracleResult]
    ],
    True,
    TestID -> "solveScenarioWithOracle: matches direct activeSetReduceSystem on 3-chain"
]
