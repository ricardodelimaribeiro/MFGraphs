(* Wolfram Language package *)
(* Disjunctive normal form reduction for systems of equations/inequalities *)
(* Solves equalities by substitution and reduces Or-branches *)

DNFReduce::usage =
"DNFReduce[xp, sys] converts the system xp && sys into disjunctive normal form
by solving equalities, substituting, and reducing branches.
DNFReduce[xp, sys, elem] is the 3-argument form that processes a single element
from a conjunction: if elem is an equality, solve and substitute into xp and sys;
otherwise, absorb elem into xp and continue with sys.";

ReplaceSolution::usage =
"ReplaceSolution[xp, sol] substitutes the Rule sol into the expression xp.
If the result has Head And, it simplifies the first conjunct.
Otherwise, it simplifies the whole expression.";

RemoveDuplicates::usage =
"RemoveDuplicates[xp] sorts (by SimplifyCount) and then DeleteDuplicates.
For example, (A && B) || (B && A) becomes (A && B) only after sorting.";

ClearSolveCache::usage =
"ClearSolveCache[] clears the internal caches used by CachedSolve and CachedReduce.
Call this between different problem instances to prevent stale results.";

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

CachedReduce[expr_, domain_] :=
  Module[{key, result},
    key = Hash[{expr, domain}, "SHA256"];
    If[KeyExistsQ[$ReduceCache, key],
      $ReduceCache[key],
      result = Reduce[expr, domain];
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

(* Conjunction: peel off the first element and process it *)
DNFReduce[xp_, And[fst_, rst_]] :=
    If[ xp === False,
        False,
        If[ Head[fst] === Or,
            (* Distribute sequentially with early exit:
               if any branch leaves xp unchanged, xp already satisfies the Or, so stop. *)
            Catch[
                Module[{accumulated = False, r},
                    Do[
                        r = DNFReduce[xp, rst, b];
                        If[r === xp, Throw[xp, "dnfAndOr"]];
                        accumulated = Or[accumulated, r],   (* Or[False, x] -> x automatically *)
                        {b, List @@ fst}
                    ];
                    If[accumulated === False, False, RemoveDuplicates[accumulated]]
                ],
                "dnfAndOr"
            ],
            (* Process the single element *)
            DNFReduce[xp, rst, fst]
        ]
    ]

(* Disjunction: reduce each branch independently, prune False results *)
DNFReduce[xp_, Or[fst_, scd_]] :=
    Module[{rfst, result1, result2},
      rfst = CachedReduce[fst, Reals];
      result1 = DNFReduce[xp, rfst];
      (* Short-circuit: if first branch already covers xp, skip second *)
      If[result1 === xp,
        result1,
        result2 = DNFReduce[xp, scd];
        Which[
          result1 === False && result2 === False, False,
          result1 === False, result2,
          result2 === False, result1,
          True, RemoveDuplicates@(Or @@ {result1, result2})
        ]
      ]
    ]

(* Fallback: absorb into the accumulated expression *)
DNFReduce[xp_, leq_] := xp && leq

(* --- DNFReduce: 3-argument forms (process one element from a conjunction) --- *)

(* Equality: solve, substitute into both xp and rst, recurse *)
DNFReduce[xp_, rst_, fst_Equal] :=
    Module[{newfst = Simplify@fst, fsol, newxp, newrst},
        If[ newfst === False,
            False,
            fsol = First@CachedSolve@newfst;
            newxp = ReplaceSolution[xp, fsol];
            If[ newxp === False,
                False,
                newrst = ReplaceSolution[rst, fsol];
                DNFReduce[newxp && fst, newrst]
            ]
        ]
    ]

(* Default: absorb the element into xp and continue with rst *)
DNFReduce[xp_, rst_, fst_] := DNFReduce[xp && fst, rst]

(* --- Backward-compatible alias --- *)
ZAnd = DNFReduce;

(* --- ReplaceSolution --- *)

ReplaceSolution[rst_?BooleanQ, sol_] := rst

ReplaceSolution[rst_, sol_] :=
    With[{newrst = rst /. sol},
        If[ Head[newrst] === And,
            And[Simplify@First@newrst, Rest@newrst],
            Simplify[newrst]
        ]
    ]

(* --- RemoveDuplicates --- *)

SortOp = SortBy[Simplify`SimplifyCount];

RemoveDuplicates[xp_And] := DeleteDuplicates[SortOp[xp]];

RemoveDuplicates[xp_Or] := DeleteDuplicates[SortOp[xp]];

RemoveDuplicates[xp_] := xp
