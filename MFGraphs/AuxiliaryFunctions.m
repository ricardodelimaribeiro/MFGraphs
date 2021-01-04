(* Wolfram Language package *)
(*TODO put usage messages outside of Private. (so we can use them) DONE*)

OtherWay::usage =
"OtherWay[{a,a\[DirectedEdge]b}] returns {b,a\[DirectedEdge]b}";

triple2path::usage =
"triple2path[{a, b, c}, G] takes 3 ordered vertices and gives the pair of edges conecting the first two and the last two. If this is not feasible, error message with {a, b, c}!";

TransitionsAt::usage = 
"";

OtherWay::usage = 
"";

AtHead::usage = 
""; 

AtTail::usage = 
"";

(*IncomingEdges::usage = 
"";

OutgoingEdges::usage = 
"";*)

(*EE::usage = 
"";

UU::usage = 
"";

RR::usage = 
"";*)
F::usage =
"for a numeric value of the current, this is the solution, m, for the H-J equation (u_x = -j)";

Intg::usage = 
"";

M::usage = 
"M[j,x] F without the Parameters association";

IntM::usage = 
"IntM[j] Intg without the Parameters association";

Begin["`Private`"]
AtHead[a_ \[DirectedEdge] b_] :=
    {b, a \[DirectedEdge] b};
(*AtHead[edge_] := {edge[[2]], edge}
AtTail[edge_] := {edge[[1]], edge}
*)
AtTail[a_ \[DirectedEdge] b_] :=
    {a, a \[DirectedEdge] b};

TransitionsAt[G_, k_] :=
    Prepend[#, k] & /@ Permutations[IncidenceList[G, k], {2}];
    

OtherWay[{c_, a_ \[DirectedEdge] b_}] :=
    {If[ c === a,
         b,
         a
     ], a \[DirectedEdge] b}
(*OtherWay[{a_, a_ \[DirectedEdge] b_}] :=
    {b, a \[DirectedEdge] b}
OtherWay[{b_, a_ \[DirectedEdge] b_}] :=
    {a, a \[DirectedEdge] b}
OtherWay[{c_, edge_}] := 
	If[c === edge[[2]], 
		AtTail[edge],
		AtHead[edge]
	]*)

triple2path[{a_, b_, c_}, G_] :=
Module[ {EL = EdgeList[G]},
    If[ SubsetQ[AdjacencyList[G, b], {a, c}],
        If[ MemberQ[EL, a \[DirectedEdge] b],
            If[ MemberQ[EL, b \[DirectedEdge] c],
                {b, a \[DirectedEdge] b, b \[DirectedEdge] c},
                {b, a \[DirectedEdge] b, c \[DirectedEdge] b}
            ],
            If[ MemberQ[EL, b \[DirectedEdge] a],
                If[ MemberQ[EL, b \[DirectedEdge] c],
                    {b, b \[DirectedEdge] a, b \[DirectedEdge] c},
                    {b, b \[DirectedEdge] a, c \[DirectedEdge] b}
                ]
            ]
        ],
        StringJoin["\nThere is no path from ", ToString[a], " to ", ToString[b], " to ", ToString[c]]
        ]
    ];

            
(*EE = Exists[#1, #2] &;

UU = Solve[#1, #2, Reals] &;

RR = Reduce[#1, #2, Reals] &;*)

F[j_?NumericQ, x_?NumericQ] :=
	First @ Values @ FindRoot[Parameters["H[x,p,m]"][x, -j/(m^(1 - Parameters["alpha"])),m], {m, 1}];
(*we are solving the h-j equation, H[x,u_x,m] == 0, for m, given that u_x = -j*m^(alpha-1) 
The hamiltonian depends on g[m]!*)


(*TODO is m becoming zero?*)
Intg[j_?NumericQ] :=
	If[j == 0.|| j == 0 , 0. ,
    -j Chop[NIntegrate[ 1/(F[j, x]^(1 - Parameters["alpha"])), {x, 0, 1}]] // Quiet
	]
(*this is the rhs of the integral of u_x from 0 to 1: ingetral_0^1 \
-j*m^(alpha-1) dx*)

M[j_?NumericQ, x_?NumericQ, edge_] :=
	First @ Values @ FindRoot[H[x, -j/(m^(1 - alpha)),m, edge], {m, 1}];


IntM[j_?NumericQ, edge_] :=	
	If[j == 0.|| j == 0 , 0. ,
    -j NIntegrate[ 1/(M[j, x, edge]^(1 - alpha)), {x, 0, 1}] // Quiet
	]


(*TODO : from previous guess, look at rhs, if m is negative, the corresponding j needs to be 0*)



(*TODO figure out what this was:
	With[{j = #}, Module [{m}, m[x_?NumericQ] := F[j, x];
		Intg[j]
	]
]&;
*)  
  
  
End[]

