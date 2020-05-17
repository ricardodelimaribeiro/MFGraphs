(* Wolfram Language package *)
MFGEquations = DataToEquations[Data];

   
jays = AssociationThread[{#,jvars[AtHead[#]] - jvars[AtTail[#]]}& /@ MFGEquations["BEL"]](*TODO check if this is what we want. or something else.*)

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
]&

world = Association[{
  
   "EqAll" ->  MFGEquations["EqAll"],
   "EqAllComp" ->  MFGEquations["EqAllComp"],
   "lhs" -> 
    AssociationThread[
   "#" -> (uvars[AtHead[#]] - uvars[AtTail[#]] -  jays[#]) & /@ MFGEquations["BEL"]], 
   "rhs" -> 
    AssociationThread[
     "#" -> (Intg(jays[#]) - jays[#]) & /@ MFGEquations["BEL"]]
   }]
 
(*edge is #.....*) 

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
		(*TODO solve and do all the things diogo said.*)
	]
]

FixedPoint[G,0]


EqValue[T_] = 
 And @@ ((T[#] == uvars[AtHead[#]] - uvars[AtTail[#]]) & /@ BEL); 
 (* TODO use the same reasoning to approach the more general cases.*)
