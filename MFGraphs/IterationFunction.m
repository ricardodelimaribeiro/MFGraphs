(* ::Package:: *)

(* Wolfram Language package *)


(*Data = DataG[2];
MFGEquations = DataToEquations[Data];*)

(*ToRules @ Or @@ Reduce[# && Data[type2] && Data[type3]]&/@Data[type1] (*if Data[type1] has head Or*)
*)
(*jays = AssociationThread[
   MFGEquations[
    "BEL"], (MFGEquations["jvars"][AtHead[#]] - 
       MFGEquations["jvars"][AtTail[#]] &) /@ MFGEquations["BEL"]];(*this is done in the variables level.*)
*)

F::usage =
"";

Intg::usage =
"";


Begin["`Private`"]
F[j_?NumericQ, x_?NumericQ] :=
    First@Values@FindRoot[Parameters["H[x,p,m]"][x, -j/(m^(1 - Parameters["alpha"])), m], {m, 1}];
(*we are solving the h-j equation for m given that u_x = -j*m^(alpha-1) *)

    
Intg[j_?NumericQ] :=
    - j NIntegrate[F[j, x]^(Parameters["alpha"] - 1), {x, 0, 1}] // Quiet
(*this is the rhs of the integral of u_x from 0 to 1: ingetral_0^1 -j*m^(alpha-1) dx*)

(*With[{j = #},
	Module [{m},
		m[x_?NumericQ] := F[j, x];
		Intg[j]
	]
]&*)

(*world = Association[{
  
   "EqAll" ->  MFGEquations["EqAll"],
   "EqAllComp" ->  MFGEquations["EqAllComp"],
   "lhs" -> 
    Association[
   (# -> MFGEquations["uvars"][AtHead[#]] - MFGEquations["uvars"][AtTail[#]] -  jays[#] &)/@ MFGEquations["BEL"]], 
   "rhs" -> 
    Association[
     (# -> Intg(jays[#]) - jays[#] &) /@ MFGEquations["BEL"]]
   }]*)
 
(*edge is #.....*) 

(*TODO build a function that returns the result for the critical congestion case*)

(*TODO learn how to document code!*)

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
(*G[world_Association] := 
Function[ jays,
	Module[{AssembledEquations},
		AssembledEquations = And@@Map[world["lhs"][#] == world["rhs"][#][jays[#]], world["edges"]];
		(*TODO solve and do all the things diogo said.*)
	]
]

FixedPoint[G,0]
*)

(*EqValue[T_] = 
 And @@ ((T[#] == uvars[AtHead[#]] - uvars[AtTail[#]]) & /@ BEL); 
*) (* TODO use the same reasoning to approach the more general cases.*)
End[]