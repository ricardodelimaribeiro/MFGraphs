#!/usr/bin/env wolframscript
(* LEGACY_GETEXAMPLEDATA: This script still depends on GetExampleData and is excluded from the active scenario-first workflow. *)
(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

(* ExactMode smoke test: verify that the new "ExactMode" -> True path on
   CriticalCongestionSolver returns determined rules in "AssoCritical" and
   the symbolic residual (when present) in "SymbolicRegion". *)

Needs["MFGraphs`"];

runCase[label_, dataExpr_, expectSymbolic_] :=
    Module[{data, d2e, res, asso, region},
        Print["--- ", label, " ---"];
        data = dataExpr;
        d2e = DataToEquations[data];
        res = CriticalCongestionSolver[d2e, "ExactMode" -> True];
        If[!AssociationQ[res],
            Print["  FAIL: result not an Association: ", res];
            Return[False, Module]
        ];
        Print["  IsFeasible      : ", IsFeasible[res]];
        Print["  Status          : ", Lookup[res, "Status", Missing["absent"]]];
        Print["  NumericBackend  : ", Lookup[res, "NumericBackendUsed", Missing["absent"]]];
        Print["  JFirstBackend   : ", Lookup[res, "JFirstBackendUsed", Missing["absent"]]];
        region = Lookup[res, "SymbolicRegion", Missing["absent"]];
        asso = Lookup[res, "AssoCritical", Missing["absent"]];
        Print["  SymbolicRegion  : ", region];
        Print["  AssoCritical    : ", asso];
        If[expectSymbolic,
            Print[If[region =!= None && !MissingQ[region], "  OK: SymbolicRegion present", "  WARN: expected SymbolicRegion"]],
            Print[If[region === None || MissingQ[region], "  OK: fully determined", "  WARN: unexpected SymbolicRegion"]]
        ];
        True
    ];

runCase["Case 12 / I1=100 (expected: fully determined)",
    GetExampleData[12] /. {I1 -> 100, U1 -> 0},
    False];
runCase["triangle with two exits (expected: symbolic region)",
    GetExampleData["triangle with two exits"] /. {I1 -> 100, U1 -> 0, S1 -> 1},
    True];

Print["Done."];

(* :!CodeAnalysis::EndBlock:: *)
