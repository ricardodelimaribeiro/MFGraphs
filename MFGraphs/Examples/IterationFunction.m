(* Wolfram Language package *)

<<<<<<< HEAD
   
jays=Association[{"edges"->jvars[AtHead[edge]] - jvars[AtTail[edge]]}]

F[j_?NumericQ, x_?NumericQ] :=
    First@Values@FindRoot[H[x, -j/(m^(1 - alpha)), m], {m, 1}];

(*TODO check if this is the right Intg function. Do the math!!*)    
    
Intg[j_] :=
    - j NIntegrate[1/(F[j, x]^(1 - alpha)), {x, 0, 1}]

(*TODO this goes in the rhs asociation. DONE*)
With[{j = #},
	Module [{m},
		m[x_?NumericQ] := F[j, x];
		Intg[j]
	]
]

world = Association[{
   "other" -> Association[EqAll && EqCompCon],
   "lhs" -> 
    Association[
     "edge" -> (uvars[AtHead[edge]] - uvars[AtTail[edge]] - 
          jvars[AtHead[edge]] - jvars[AtTail[edge]]) & /@ 
      EdgeList[BG]], 
   "rhs" -> 
    Association[
     "edge" -> (Intg(jays[edge])) & /@ EdgeList[BG]]
   }]

(*TODO build a function that returns the result for the critical congestion case*)


(*TODO make the iterative function work on all the values.*)

(*G[world_Association] := 
Function[ jays,  (*jays associates numeric values to the edges \
of BG. Run "the code" in the critical congestion and feed it to G.*)
    Module[{ItValues,EqItAll, EqAllRNZ, ResolvedRNZ, URNZ, JTJNZ, jvars, uvars, BG, js, jts, us, EqNoInt},
    ItValues = 
        And @@ ((Intg[jays[AtHead[#]] - jays[AtTail[#]]] == 
        uvars[AtHead[#]] - uvars[AtTail[#]]) & /@ EdgeList[BG]);
    EqItAll = EqNoInt && ItValues;
    Print["Reducing ..."];
    EqAllRNZ = Reduce[EqItAll, Reals];
    Print["Eliminating js and jts ..."];
    ResolvedRNZ = Resolve[EE[Join[js, jts], EqAllRNZ], Reals];
    Print["Finding value function rules ... "];
    URNZ = UU[ResolvedRNZ, us] // Flatten;
    Print["Finding current and transition current rules ..."];
    JTJNZ = UU[EqAllRNZ /. URNZ, Join[jts, js]] // Flatten;
    jvars /. JTJNZ
    ]
    (*TODO what do we do when there are multiple solutions*)
]
DataToEquations[Data];
*)

(*TODO finish this which diogo started.*)
G[world_Association] := 
Function[ jays,
	Module[{AssembledEquations},
		AssembledEquations = And@@Map[world["lhs"][#] == world["rhs"][#][jays[#]], world["edges"]];
	]
]

FixedPoint[G,0]


EqValue[T_] = 
 And @@ ((T[#] == uvars[AtHead[#]] - uvars[AtTail[#]]) & /@ BEL); 
 (* TODO use the same reasoning to approach the more general cases.*)
=======
F[j_?NumericQ, x_?NumericQ] :=
    First@Values@FindRoot[H[x, -j/(m^(1 - alpha)), m], {m, 1}];

(*TODO check if this is the right Intg function. Do the math!!*)    
    
Intg[j_] :=
    - j NIntegrate[1/(F[j, x]^(1 - alpha)), {x, 0, 1}]

(*TODO this goes in the rhs asociation*)
With[{j = #},
	Module [{m},
		m[x_?NumericQ] := F[j, x];
		Intg[j]
	]
]



(*TODO build a function that returns the result for the critical congestion case*)


(*TODO make the iterative function work on all the values.*)

G[world_Association] := 
Function[ jays,  (*jays associates numeric values to the edges \
of BG. Run "the code" in the critical congestion and feed it to G.*)
    Module[{ItValues,EqItAll, EqAllRNZ, ResolvedRNZ, URNZ, JTJNZ, jvars, uvars, BG, js, jts, us, EqNoInt},
    ItValues = 
        And @@ ((Intg[jays[AtHead[#]] - jays[AtTail[#]]] == 
        uvars[AtHead[#]] - uvars[AtTail[#]]) & /@ EdgeList[BG]);
    EqItAll = EqNoInt && ItValues;
    Print["Reducing ..."];
    EqAllRNZ = Reduce[EqItAll, Reals];
    Print["Eliminating js and jts ..."];
    ResolvedRNZ = Resolve[EE[Join[js, jts], EqAllRNZ], Reals];
    Print["Finding value function rules ... "];
    URNZ = UU[ResolvedRNZ, us] // Flatten;
    Print["Finding current and transition current rules ..."];
    JTJNZ = UU[EqAllRNZ /. URNZ, Join[jts, js]] // Flatten;
    jvars /. JTJNZ
    ](*TODO what do we do when there are multiple solutions*)
]
DataToEquations[Data];



(*TODO finish this which diogo started.*)
G[world_Association] := 
Function[ jays,
	Module[{AssembledEquations},
		AssembledEquations = And@@Map[world["LHS"][#] == world["rhs"][#][jays[#]], world["edges"]];
	]
]

FixedPoint[G,0]


EqValue[T_] = 
 And @@ ((T[#] == uvars[AtHead[#]] - uvars[AtTail[#]]) & /@ BEL); (* TODO use the same reasoning to approach the more general cases.*)
>>>>>>> branch 'master' of https://github.com/ricardodelimaribeiro/MFGraphs
(*assemble all the equations out of Eq function! then define Eq*)

(*TODO Write the current values as an association. then fix this: Reduce[EqNoInt&&EqValue[], Reals] // AbsoluteTiming*)

(*TODO decide how to iterate the system. Make everything into associations!*)

