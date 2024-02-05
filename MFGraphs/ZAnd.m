(* Wolfram Language package *)

ZAnd::usage =
"
ZAnd[xp,xps] returns a system which is equivalent to xp&&xps in disjunctive normal form. 
ZAnd[xp, And[fst, scd] 
ZAnd[xp, eq] returns xp with the solution of eq replaced in it together with eq.
ZAnd[xp, Or[fst,scd]] returns ZAnd[xp, fst]||ZAnd[xp, scd]";

ZAnd[_, False] = 
	False

ZAnd[False, _] = 
	False

ZAnd[xp_, True] := 
	xp

ZAnd[xp_, eq_Equal] := 
 With[{sol = Solve[eq]}, 
  If[sol === {}, 
   False, (xp /. First@sol) && And @@ (First@sol /. Rule -> Equal) // 
    Simplify]]

ZAnd[xp_, orxp_Or] := (ZAnd[xp, #] & /@ orxp) // RemoveDuplicates

ZAnd[xp_, leq_] := Simplify[xp && leq]

ZAnd[xp_, andxp_And] := 
 With[{fst = First[andxp], rst = Rest[andxp]}, 
  Print["ZAnd: head is and: ", fst];
  Which[
  	Head[fst] === Or, 
  		ReZAnd[xp, rst] /@ fst // RemoveDuplicates, 
  	True, 
   		ReZAnd[xp, rst, fst]
  ]
 ]


(*Operator form of ReZAnd*)
ReZAnd[xp_, rst_] := ReZAnd[xp, rst, #] &

ReZAnd[xp_, rst_, fst_Equal] := 
 Module[{fsol = First@Solve@fst // Quiet, newrst, newxp}, 
  newrst = ReplaceSolution[rst, fsol];
  newxp = xp /. fsol;
  Print["ReZAnd: ", fst, " ", newrst, " ", newxp];
  ZAnd[newxp && fst, newrst]]

ReZAnd[xp_, rst_, fst_] := ZAnd[xp && fst, rst]

ReplaceSolution::usage =
"ReplaceSolution[xp,sol] substitutes the Rule, solution, on the expression xp.
After that, it reduces each sub-expression over the Reals if the Head of the expression is And. 
Otherwise, it reduces the whole expression over the Reals.";

ReplaceSolution[rst_?BooleanQ, sol_] := rst

ReplaceSolution[rst_, sol_] := Module[{newrst}, newrst = rst /. sol;
  If[Head[newrst] === And, Reduce[#, Reals] & /@ newrst, 
   Reduce[newrst, Reals]]]

(*ReplaceSolution[rst_,sol_]:=Simplify[rst/. sol]*)
RemoveDuplicates::usage =
"RemoveDuplicates[xp] sorts and then DeleteDuplicates. 
We need to sort because DeleteDuplicates only deletes identical expressions.
For example (A&&B)||(B&&A) becomes (A&&B) only after sorting.
"
RemoveDuplicates[xp_And] := DeleteDuplicates[Sort[xp]];
RemoveDuplicates[xp_Or] := DeleteDuplicates[Sort[xp]];
RemoveDuplicates[xp_] := xp
