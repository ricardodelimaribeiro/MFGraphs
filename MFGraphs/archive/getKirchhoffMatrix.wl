(* Wolfram Language package *)
(*
   getKirchhoffMatrix.wl — archived from systemTools.wl.

   Wrapped getKirchhoffLinearSystem with a placeholder cost-function slot
   intended for future congestion-cost injection. Zero in-tree callers at
   the time of archival; the placeholder body returned a constant Function
   producing {0&, 0&, ...} regardless of input mass.

   Restoration: re-add the ::usage declaration and definition to
   systemTools.wl alongside getKirchhoffLinearSystem.
*)

BeginPackage["getKirchhoffMatrixArchive`", {"systemTools`"}]

getKirchhoffMatrix::usage =
"getKirchhoffMatrix[sys] returns the entry current vector, Kirchhoff matrix, \
(critical congestion) cost function placeholder, and the variables in the \
order corresponding to the Kirchhoff matrix. The third slot is retained for \
backward compatibility.";

Begin["`Private`"]

getKirchhoffMatrix[sys_] :=
    Module[{b, km, vars, cost, cCost},
        {b, km, vars} = getKirchhoffLinearSystem[sys];
        If[Length[vars] === 0,
            Return[{b, km, <||>, vars}, Module]
        ];
        cost  = AssociationThread[vars, (0&) /@ vars];
        cCost = cost /@ vars /. MapThread[Rule, {vars, #}]&;
        {b, km, cCost, vars}
    ];

End[]

EndPackage[]
