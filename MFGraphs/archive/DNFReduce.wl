(* Wolfram Language package *)
(*
   DNFReduce: Symbolic logic reduction for MFG systems.

   This module handles the combinatorial complexity of network Mean Field Games
   by converting constraints into Disjunctive Normal Form (DNF). it simplifies
   systems by solving equalities via substitution and pruning redundant
   logical branches through subsumption checking.
*)

BeginPackage["MFGraphs`"];

DNFReduce::usage =
"DNFReduce[xp, sys] converts the system xp && sys into disjunctive normal form
by solving equalities, substituting, and reducing branches.
DNFReduce[xp, sys, elem] is the 3-argument form that processes a single element
from a conjunction: if elem is an equality, solve and substitute into xp and sys;
otherwise, absorb elem into xp and continue with sys.";

SubstituteSolution::usage =
"SubstituteSolution[xp, sol] substitutes the Rule sol into the expression xp.
If the result has Head And, it simplifies the first conjunct.
Otherwise, it simplifies the whole expression.";

DeduplicateByComplexity::usage =
"DeduplicateByComplexity[xp] sorts (by SimplifyCount) and then DeleteDuplicates.
For example, (A && B) || (B && A) becomes (A && B) only after sorting.";

ClearSolveCache::usage =
"ClearSolveCache[] clears the internal caches used by CachedSolve and CachedReduce.
Call this between different problem instances to prevent stale results.";

ReduceDisjuncts::usage =
"ReduceDisjuncts[expr] reduces the number of disjuncts in a DNF expression by
removing subsumed branches. For small expressions (<= $ReduceDisjunctsThreshold
disjuncts) uses Reduce directly; for larger expressions uses subsumption pruning
followed by Reduce on the survivors.";

$ReduceDisjunctsThreshold::usage =
"$ReduceDisjunctsThreshold controls the cutoff between Full Reduce and Subsumption
strategies in ReduceDisjuncts. Default is 200.";
$ReduceDisjunctsThreshold = 300;

RemoveDuplicates::usage = "RemoveDuplicates is a backward-compatibility alias for DeduplicateByComplexity.";
ReplaceSolution::usage = "ReplaceSolution is a backward-compatibility alias for SubstituteSolution.";

Begin["`Private`"];

(* --- Memoization caches --- *)

$SolveCache = <||>;
$ReduceCache = <||>;

CachedSolve[eq_] :=
  Module[{key, result},
    key = Hash[eq, "SHA256"];
    If[KeyExistsQ[$SolveCache, key],
      $SolveCache[key],
      result = Solve[eq];
      $SolveCache[key] = result;
      result
    ]
  ];

CachedReduce[expr_, vars_, domain_] :=
  Module[{key, result},
    key = Hash[{expr, vars, domain}, "SHA256"];
    If[KeyExistsQ[$ReduceCache, key],
      $ReduceCache[key],
      result = Reduce[expr, vars, domain];
      $ReduceCache[key] = result;
      result
    ]
  ];

ClearSolveCache[] := ($SolveCache = <||>; $ReduceCache = <||>);

(* --- DNFReduce: 2-argument forms (main dispatch) --- *)

DNFReduce[_, False] := False

DNFReduce[False, _] := False

DNFReduce[xp_, True] := xp

(* Single equality: delegate to 3-arg form *)
DNFReduce[xp_, eq_Equal] := DNFReduce[xp, True, eq]

(* Helper: flatten post-substitution results back to a list of conjuncts *)
FlattenConjuncts[True] := {}
FlattenConjuncts[False] := {False}
FlattenConjuncts[x_And] := List @@ x
FlattenConjuncts[x_] := {x}

(* Conjunction: flatten to list and process iteratively to avoid deep recursion.
   The original recursive version peeled off one conjunct at a time, causing stack
   overflow on large networks like Grid1020 (200 vertices). This iterative version
   processes conjuncts in a While loop, only recursing for Or-distribution. *)
DNFReduce[xp_, sys_And] :=
    Module[{conjuncts = List @@ sys, xpAcc = xp, i = 1, elem,
            newfst, sol, fsol, newxp, rest, newRest},
        While[i <= Length[conjuncts],
            If[xpAcc === False, Return[False, Module]];
            elem = conjuncts[[i]];

            Which[
                elem === True,
                    i++,

                elem === False,
                    Return[False, Module],

                Head[elem] === Or,
                    (* Distribute: process each Or-branch with remaining conjuncts.
                       This is the only place we recurse, via the existing 3-arg forms. *)
                    rest = If[i >= Length[conjuncts], True, And @@ conjuncts[[i + 1 ;;]]];
                    Return[
                        Catch[
                            Module[{accumulated = False, r},
                                Do[
                                    r = DNFReduce[xpAcc, rest, b];
                                    If[r === xpAcc, Throw[xpAcc, "dnfAndOr"]];
                                    accumulated = Or[accumulated, r],
                                    {b, List @@ elem}
                                ];
                                If[accumulated === False, False, DeduplicateByComplexity[accumulated]]
                            ],
                            "dnfAndOr"
                        ],
                        Module
                    ],

                Head[elem] === Equal,
                    (* Equality: solve, substitute into xpAcc and all remaining conjuncts *)
                    newfst = Simplify @ elem;
                    If[newfst === False, Return[False, Module]];
                    sol = Quiet[CachedSolve @ newfst];
                    fsol = If[MatchQ[sol, {{__Rule}, ___}], First[sol], {}];
                    newxp = SubstituteSolution[xpAcc, fsol];
                    If[newxp === False, Return[False, Module]];
                    xpAcc = newxp && elem;
                    (* Substitute into all remaining conjuncts *)
                    If[i < Length[conjuncts],
                        newRest = SubstituteSolution[And @@ conjuncts[[i + 1 ;;]], fsol];
                        conjuncts = Join[conjuncts[[;; i]], FlattenConjuncts[newRest]];
                    ];
                    i++,

                True,
                    (* Default: absorb into xpAcc *)
                    xpAcc = xpAcc && elem;
                    i++
            ]
        ];
        xpAcc
    ]

(* Disjunction: reduce each branch independently, prune False results *)
DNFReduce[xp_, Or[fst_, scd_]] :=
    Module[{rfst, result1, result2},
      rfst = CachedReduce[fst, Variables[fst], Reals];
      result1 = DNFReduce[xp, rfst];
      (* Short-circuit: if first branch already covers xp, skip second *)
      If[result1 === xp,
        result1,
        result2 = DNFReduce[xp, scd];
        Which[
          result1 === False && result2 === False, False,
          result1 === False, result2,
          result2 === False, result1,
          True, DeduplicateByComplexity@(Or @@ {result1, result2})
        ]
      ]
    ]

(* Fallback: absorb into the accumulated expression *)
DNFReduce[xp_, leq_] := xp && leq

(* --- DNFReduce: 3-argument forms (process one element from a conjunction) --- *)

(* Equality: solve, substitute into both xp and rst, recurse *)
DNFReduce[xp_, rst_, fst_Equal] :=
    Module[{newfst = Simplify@fst, fsol, sol, newxp, newrst},
        If[ newfst === False,
            False,
            sol = Quiet[CachedSolve@newfst];
            fsol = If[MatchQ[sol, {{__Rule}, ___}], First[sol], {}];
            newxp = SubstituteSolution[xp, fsol];
            If[ newxp === False,
                False,
                newrst = SubstituteSolution[rst, fsol];
                DNFReduce[newxp && fst, newrst]
            ]
        ]
    ]

(* Default: absorb the element into xp and continue with rst *)
DNFReduce[xp_, rst_, fst_] := DNFReduce[xp && fst, rst]

(* --- SubstituteSolution --- *)

SubstituteSolution[rst_?BooleanQ, sol_] := rst

SubstituteSolution[rst_, sol_] :=
    With[{newrst = rst /. sol},
        If[ Head[newrst] === And,
            And[Simplify@First@newrst, Rest@newrst],
            Simplify[newrst]
        ]
    ]

(* --- DeduplicateByComplexity --- *)
(* Sort operands of And/Or for canonical form, then remove duplicates.
   Sort uses Mathematica's canonical ordering, which is faster than
   SortBy[SimplifyCount] and provides a unique canonical form. *)

DeduplicateByComplexity[xp_] := DeleteDuplicates[Sort[xp]] /; MatchQ[Head[xp], And|Or]

DeduplicateByComplexity[xp_] := xp

(* --- Subsumption pruning --- *)
(* Remove branch j if some branch i implies it (i.e., i is more general). *)

SubsumptionPrune[branches_List, perCheckTimeout_:2] :=
  Module[{n = Length[branches], keep, result},
    keep = ConstantArray[True, n];
    Do[
      If[!keep[[i]], Continue[]];
      Do[
        If[!keep[[j]], Continue[]];
        (* Does branch i imply branch j? If so, j is redundant. *)
        Module[{impl = Implies[branches[[i]], branches[[j]]]},
          result = Quiet@TimeConstrained[
            Reduce[impl, Variables[impl], Reals],
            perCheckTimeout, $TimedOut
          ]
        ];
        If[result === True,
          keep[[j]] = False;
          Continue[]
        ];
        (* Does branch j imply branch i? If so, i is redundant. *)
        Module[{impl = Implies[branches[[j]], branches[[i]]]},
          result = Quiet@TimeConstrained[
            Reduce[impl, Variables[impl], Reals],
            perCheckTimeout, $TimedOut
          ]
        ];
        If[result === True,
          keep[[i]] = False;
          Break[]
        ],
        {j, i + 1, n}
      ],
      {i, 1, n}
    ];
    Pick[branches, keep]
  ];

(* --- ReduceDisjuncts: tiered post-reduction of DNF output --- *)

ReduceDisjuncts[expr_] := expr /; Head[expr] =!= Or

ReduceDisjuncts[expr_Or] :=
  Module[{branches, n, reduced, time},
    branches = List @@ expr;
    n = Length[branches];
    MFGPrint["ReduceDisjuncts: ", n, " disjuncts, ", LeafCount[expr], " leaves"];
    If[n <= $ReduceDisjunctsThreshold,
      (* Small: Full Reduce can handle it directly *)
      {time, reduced} = AbsoluteTiming[
        Quiet@TimeConstrained[Reduce[expr, Variables[expr], Reals], 120, $TimedOut]
      ];
      If[reduced === $TimedOut,
        MFGPrint["ReduceDisjuncts: Full Reduce timed out, falling back to subsumption"];
        reduced = SubsumptionPrune[branches];
        reduced = Or @@ reduced,
        MFGPrint["ReduceDisjuncts: Full Reduce ", n, " -> ",
          If[Head[reduced] === Or, Length[reduced], 1], " in ", time, "s"]
      ],
      (* Large: subsumption first, then Reduce on survivors *)
      {time, reduced} = AbsoluteTiming[SubsumptionPrune[branches]];
      MFGPrint["ReduceDisjuncts: Subsumption ", n, " -> ", Length[reduced], " in ", time, "s"];
      If[Length[reduced] > 1 && Length[reduced] <= $ReduceDisjunctsThreshold,
        (* Survivors are small enough for Reduce *)
        Module[{r2, t2},
          {t2, r2} = AbsoluteTiming[
            Quiet@TimeConstrained[Reduce[Or @@ reduced, Variables[Or @@ reduced], Reals], 120, $TimedOut]
          ];
          If[r2 =!= $TimedOut,
            MFGPrint["ReduceDisjuncts: Reduce on survivors -> ",
              If[Head[r2] === Or, Length[r2], 1], " in ", t2, "s"];
            reduced = r2,
            reduced = Or @@ reduced
          ]
        ],
        reduced = Or @@ reduced
      ]
    ];
    reduced
  ];

(* --- Backward compatibility aliases --- *)
RemoveDuplicates = DeduplicateByComplexity;
ReplaceSolution = SubstituteSolution;

End[];

EndPackage[];
