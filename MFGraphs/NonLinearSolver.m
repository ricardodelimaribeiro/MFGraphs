(*Wolfram Language package*)
Get["/Users/ribeirrd/eclipse-workspace/MFGraphs/MFGraphs/D2E2.m"];
Clear[H,Cost];
NonLinear::usage =
    "NonLinear[Eqs] takes an association resulting from Data2Equations and returns an approximation to the solution of the non-critical congestion case with alpha = value, specified by alpha[edge_]:= value.";
IsNonLinearSolution::usage = 
"IsNonLinearSolution[Eqs] extracts AssoNonCritial and checks equations and inequalities. The right and left hand sides of the nonlinear equations are shown with the sup-norm of the difference."
NonLinear[Eqs_] :=
    Module[ {AssoCritical, PreEqs = Eqs, AssoNonCritical, NonCriticalList, js = 1, MaxIter = 15},
        If[ KeyExistsQ[PreEqs, "AssoCritical"], 
            (*if there is already an approximation for the non congestion case, use it!*)
            AssoNonCritical = Lookup[PreEqs, "AssoNonCritical", PreEqs["AssoCritical"]],
            PreEqs = MFGPreprocessing[PreEqs];
            js = Lookup[PreEqs, "js", $Failed];
            AssoNonCritical = AssociationThread[js, 0 js]
        ];
        NonCriticalList = FixedPointList[NonLinearStep[PreEqs], AssoNonCritical, MaxIter];
        Print["Iterated ", Length[NonCriticalList]-1, " times out of ", MaxIter];
        AssoCritical = Lookup[PreEqs, "AssoCritical", NonCriticalList[[2]]];
        AssoNonCritical = NonCriticalList // Last;
        Join[PreEqs, Association[{"AssoCritical" -> AssoCritical, "AssoNonCritical" -> AssoNonCritical}]]
    ];
style = Style[#,Bold,Blue]&;
NonLinearStep[Eqs_][approxSol_] :=
    Module[ {approxJs, approx, js, Nrhs, Nlhs, Newlhs, Newrhs},
        js = Lookup[Eqs, "js", $Failed];
        Nrhs = Lookup[Eqs, "Nrhs", $Failed];
        Nlhs = Lookup[Eqs, "Nlhs", $Failed];
        approxJs = KeyTake[approxSol, js];
        approx = MFGSystemSolver[Eqs][approxJs];
        Newlhs = N[Nlhs/.approx];
        Newrhs = Nrhs/.approx;
        Print[Newlhs,"\n",Newrhs];
        Print[style@"Max error for non-linear solution: ", Norm[Newlhs-Newrhs,Infinity]];
        approx
    ];
    
IsNonLinearSolution[Eqs_] :=
    Module[ {EqEntryIn, EqValueAuxiliaryEdges, EqSwitchingByVertex, EqCompCon, 
        EqBalanceSplittingCurrents, EqCurrentCompCon, EqTransitionCompCon,
        EqPosJs, EqPosJts, Nrhs, Nlhs, style, styler, styleg, assoc},
        assoc = Lookup[Eqs, "AssoNonCritical", <||>];
        (*Print[assoc];*)
        style = Style[#,Bold,Blue]&;
        styler = Style[#,Bold,Red]&;
        styleg = Style[#,Bold,Darker@Green]&;
        EqEntryIn = And @@ Lookup[Eqs, "EqEntryIn", $Failed];
        EqCompCon = Lookup[Eqs, "EqCompCon", $Failed];
        EqCurrentCompCon = Lookup[Eqs, "EqCurrentCompCon", $Failed];
        EqTransitionCompCon = Lookup[Eqs, "EqTransitionCompCon", $Failed];
        EqPosJs = Lookup[Eqs, "EqPosJs", $Failed];
        EqPosJts = Lookup[Eqs, "EqPosJts", $Failed];
        EqSwitchingByVertex = And @@ Lookup[Eqs, "EqSwitchingByVertex", $Failed];
        EqBalanceSplittingCurrents = Lookup[Eqs, "EqBalanceSplittingCurrents", $Failed];
        EqValueAuxiliaryEdges = Lookup[Eqs, "EqValueAuxiliaryEdges", $Failed];
        Nrhs = Lookup[Eqs, "Nrhs", $Failed];
        Nlhs = Lookup[Eqs, "Nlhs", $Failed];
        bool = (EqEntryIn&&EqValueAuxiliaryEdges&&EqCompCon&&
            EqBalanceSplittingCurrents&&EqCurrentCompCon&&EqTransitionCompCon&&
            EqPosJs&&EqPosJts&&EqSwitchingByVertex)/.assoc;
        If[ bool,
            Print["All restrictions are ", styleg@"True"],
            Print["At least one of the restrictions is ", styler@bool];
            Print[style@"EqEntryIn: ", EqEntryIn/.assoc, "\n", EqEntryIn];
            Print[style@"EqValueAuxiliaryEdges: ", EqValueAuxiliaryEdges/.assoc, "\n", EqValueAuxiliaryEdges];
            Print[style@"EqCompCon: ", EqCompCon/.assoc, "\n", EqCompCon];
            Print[style@"EqBalanceSplittingCurrents: ", EqBalanceSplittingCurrents/.assoc, "\n", EqBalanceSplittingCurrents];
            Print[style@"EqCurrentCompCon: ", EqCurrentCompCon/.assoc, "\n",EqCurrentCompCon];
            Print[style@"EqTransitionCompCon: ", EqTransitionCompCon/.assoc, "\n", EqTransitionCompCon];
            Print[style@"EqPosJs: ", EqPosJs/.assoc, "\n", EqPosJs];
            Print[style@"EqPosJts: ", EqPosJts/.assoc, "\n", EqPosJts];
            Print[style@"EqSwitchingByVertex: ", EqSwitchingByVertex/.assoc, "\n", EqSwitchingByVertex]
        ];
        Print[style@"Nlhs: ", N[Nlhs/.assoc], "\n", Nlhs];
        Print[style@"Nrhs: ", N[Nrhs/.assoc], "\n", Nrhs];
        Print[style@"Max error for non-linear solution: ", Norm[N[Nlhs/.assoc]-(Nrhs/.assoc),Infinity]];
        N@assoc
    ];
    
RoundValues[x_?NumberQ] :=
    Round[x, 10^-10]

RoundValues[Rule[a_, b_]] :=
    Rule[a, RoundValues[b]]

RoundValues[x_List] :=
    RoundValues /@ x

RoundValues[x_Association] :=
    RoundValues /@ x

alpha[edge_] :=
    1

g[m_, edge_] :=
    m

A = 0.5;

V::usage = 
""
V = Function[{x, edge}, W[x, A]];

W::usage =
""
W = Function[{y, a}, a Sin[2 Pi (y + 1/4)]^2];

H::usage = (*TODO include the other parameter here, if we call it beta, change the the other beta!*)
"H[xi,p,m, edge] is the Hamiltonian function for the edges.
edge is an edge from the graph, i.e. DirectedEdge[r,i]."
(*H = Function[{xi,p,m, edge}, p^2/(2 m^alpha[edge]) + V[xi, edge] - g[m]];*)
H = Function[{xi,p,m, edge}, p^2/(2 m^alpha[edge]) + V[xi, edge] - g[m, edge]];

(*H[xi_,p_,m_, edge_]:= p^2/(2 m^alpha[edge]) + V[xi, edge] - g[m, edge];*)
U::usage =
"U[x, edge, Eqs, sol]"
U[x_?NumericQ , edge_, Eqs_Association, sol_] :=
    Module[ {jay, uT},
        jay = Eqs["SignedCurrents"][edge]/.sol;
        uT = Eqs["uvars"][AtTail[edge]]/.sol;
        If[ jay == 0.|| jay == 0,
            uT,
            uT - jay NIntegrate[M[jay, y, edge]^(alpha[edge] - 1), {y, 0, x}]
        ]
    ]

M[j_?NumericQ, x_?NumericQ, edge_] :=
    If[ j == 0. || j == 0,
        0.,
        Values@First@FindRoot[H[x, -j m^(alpha[edge] - 1), m, edge], {m, 1}]
    ];

IntM[j_?NumericQ, edge_] :=
    If[ j == 0. || j == 0,
        0,
        j NIntegrate[M[j, x, edge]^(alpha[edge] - 1), {x, 0, 1}] // Quiet
    ];

Cost[j_, edge_] :=
    IntM[j, edge];

plotM[Eqs_, string_, edge_] :=
    Module[ {jays, sol},
        sol = Eqs[string];
        jays = Eqs["SignedCurrents"][edge] /. sol;
        sol = Eqs[string];
        Plot[M[jays /. sol, x, edge], {x, 0, 1}, 
         PlotLabel -> edge,(*PlotRange\[Rule]{-0.1,2.4},*)
         GridLines -> Automatic]
    ]
   
plotMs[Eqs_, string_] :=
    plotM[Eqs, string, #] & /@ crit["BEL"];
    
plotU[Eqs_, string_, edge_] :=
    Module[ {sol},
        sol = Eqs[string];
        Plot[U[x, edge, Eqs, sol], {x, 0, 1}, PlotLabel -> edge(*,
         PlotRange\[Rule]{1-0.1,4.7}*), GridLines -> Automatic]
    ];

plotUs[Eqs_, string_] :=
    plotU[Eqs, string, #] & /@ Eqs["BEL"]