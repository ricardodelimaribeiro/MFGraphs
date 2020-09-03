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
"";

Intg::usage = 
"";

Begin["`Private`"]
AtHead[a_ \[DirectedEdge] b_] :=
    {b, a \[DirectedEdge] b};

AtTail[a_ \[DirectedEdge] b_] :=
    {a, a \[DirectedEdge] b};

TransitionsAt[G_, k_] :=
    Prepend[#, k] & /@ Permutations[IncidenceList[G, k], {2}];
    

OtherWay[{c_, a_ \[DirectedEdge] b_}] :=
    {If[ c === a,
         b,
         a
     ], a \[DirectedEdge] b}


triple2path[{a_, b_, c_}, G_] :=
    Module[ {EL, FG = G},
        EL = EdgeList[FG];
        If[ SubsetQ[AdjacencyList[G, b], {a, c}],
            If[ MemberQ[EL, a \[DirectedEdge] b],
                If[ MemberQ[EL, b \[DirectedEdge] c],
                    {b, a \[DirectedEdge] b, b \[DirectedEdge] c},
                    {b, a \[DirectedEdge] b, 
                    c \[DirectedEdge] b}
                ],
                If[ MemberQ[EL, b \[DirectedEdge] a],
                    If[ MemberQ[EL, b \[DirectedEdge] c],
                        {b, b \[DirectedEdge] a, 
                        b \[DirectedEdge] c},
                        {b, b \[DirectedEdge] a, 
                        c \[DirectedEdge] b}
                    ]
                ]
            ],
            StringJoin["\n There is no path from ", ToString[a], " to ", 
             ToString[b], " to ", ToString[c]]
        ]
    ];

            
(*EE = Exists[#1, #2] &;

UU = Solve[#1, #2, Reals] &;

RR = Reduce[#1, #2, Reals] &;*)

F[j_?NumericQ, x_?NumericQ] :=
	First @ Values @ FindRoot[Parameters["H[x,p,m]"][x, -j/(m^(1 - Parameters["alpha"])),m], {m, 1}];
(*we are solving the h-j equation for m given that u_x = -j*m^(alpha-1) *)


(*TODO is m becoming zero?*)
Intg[j_?NumericQ] :=
	If[j == 0. , 0. ,
    - j NIntegrate[ 1/(F[j, x]^(1 - Parameters["alpha"])), {x, 0, 1}] // Quiet
	]
(*this is the rhs of the integral of u_x from 0 to 1: ingetral_0^1 \
-j*m^(alpha-1) dx*)

(*TODO : from previous guess, look at rhs, if m is negative, the corresponding j needs to be 0*)



(*TODO figure out what this was:
	With[{j = #}, Module [{m}, m[x_?NumericQ] := F[j, x];
		Intg[j]
	]
]&;
*)  
  
  
End[]

