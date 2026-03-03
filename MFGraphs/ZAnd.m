(* Wolfram Language package *)
(* Boolean algebra utilities for disjunctive normal form reduction *)
(* Canonical source for ZAnd, ReZAnd, ReplaceSolution, RemoveDuplicates *)

ZAnd::usage =
"ZAnd[xp, xps] returns a system equivalent to xp && xps in disjunctive normal form.
ZAnd[xp, And[fst, scd]] processes the conjunction by handling Or branches via ReZAnd.
ZAnd[xp, eq] returns xp with the solution of eq replaced in it, together with eq.
ZAnd[xp, Or[fst, scd]] returns ZAnd[xp, fst] || ZAnd[xp, scd].";

ReZAnd::usage =
"ReZAnd[xp, rst, fst] is the recursive helper for ZAnd.
ReZAnd[xp, rst] returns the operator form ReZAnd[xp, rst, #]&.";

ReplaceSolution::usage =
"ReplaceSolution[xp, sol] substitutes the Rule sol into the expression xp.
If the result has Head And, it simplifies the first conjunct.
Otherwise, it simplifies the whole expression.";

RemoveDuplicates::usage =
"RemoveDuplicates[xp] sorts (by SimplifyCount) and then DeleteDuplicates.
Sorting is needed because DeleteDuplicates only removes identical expressions.
For example, (A && B) || (B && A) becomes (A && B) only after sorting.";

ClearSolveCache::usage =
"ClearSolveCache[] clears the internal cache used by CachedSolve.
Call this between different problem instances to prevent stale results.";

(* --- Solve memoization cache --- *)

$SolveCache = <||>;

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

ClearSolveCache[] := ($SolveCache = <||>);

(* --- ZAnd: disjunctive normal form converter --- *)

ZAnd[_, False] := False

ZAnd[False, _] := False

ZAnd[xp_, True] := xp

ZAnd[xp_, eq_Equal] := ReZAnd[xp, True, eq]

ZAnd[xp_, And[fst_, rst_]] :=
    If[ xp === False,
        False,
        If[ Head[fst] === Or,
            RemoveDuplicates@(ReZAnd[xp, rst] /@ fst),
            ReZAnd[xp, rst, fst]
        ]
    ]

ZAnd[xp_, Or[fst_, scd_]] :=
    Module[{rfst, result1, result2},
      rfst = Reduce[fst, Reals];
      result1 = ZAnd[xp, rfst];
      result2 = ZAnd[xp, scd];
      (* Skip False branches in the Or result *)
      Which[
        result1 === False && result2 === False, False,
        result1 === False, result2,
        result2 === False, result1,
        True, RemoveDuplicates@(Or @@ {result1, result2})
      ]
    ]

ZAnd[xp_, leq_] := xp && leq

(* --- ReZAnd: recursive helper --- *)

(* Operator form *)
ReZAnd[xp_, rst_] := ReZAnd[xp, rst, #] &

(* Equality case: solve and substitute, with memoized Solve *)
ReZAnd[xp_, rst_, fst_Equal] :=
    Module[{newfst = Simplify@fst},
        If[ newfst === False,
            False,
            With[{fsol = First@CachedSolve@newfst},
                ZAnd[ReplaceSolution[xp, fsol] && fst, ReplaceSolution[rst, fsol]]
            ]
        ]
    ]

(* Default case: absorb into the first argument *)
ReZAnd[xp_, rst_, fst_] := ZAnd[xp && fst, rst]

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
