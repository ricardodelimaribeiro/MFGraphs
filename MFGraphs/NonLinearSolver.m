(*Wolfram Language package*)
(*Get["/Users/ribeirrd/eclipse-workspace/MFGraphs/MFGraphs/D2E2.m"];*)
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
    Module[ {EqEntryIn, EqValueAuxiliaryEdges, IneqSwitchingByVertex, AltOptCond, 
        EqBalanceSplittingFlows, AltFlows, AltTransitionFlows,
        IneqJs, IneqJts, Nrhs, Nlhs, style, styler, styleg, assoc},
        assoc = Lookup[Eqs, "AssoNonCritical", <||>];
        (*Print[assoc];*)
        style = Style[#,Bold,Blue]&;
        styler = Style[#,Bold,Red]&;
        styleg = Style[#,Bold,Darker@Green]&;
        EqEntryIn = And @@ Lookup[Eqs, "EqEntryIn", $Failed];
        AltOptCond = Lookup[Eqs, "AltOptCond", $Failed];
        AltFlows = Lookup[Eqs, "AltFlows", $Failed];
        AltTransitionFlows = Lookup[Eqs, "AltTransitionFlows", $Failed];
        IneqJs = Lookup[Eqs, "IneqJs", $Failed];
        IneqJts = Lookup[Eqs, "IneqJts", $Failed];
        IneqSwitchingByVertex = And @@ Lookup[Eqs, "IneqSwitchingByVertex", $Failed];
        EqBalanceSplittingFlows = Lookup[Eqs, "EqBalanceSplittingFlows", $Failed];
        EqValueAuxiliaryEdges = Lookup[Eqs, "EqValueAuxiliaryEdges", $Failed];
        Nrhs = Lookup[Eqs, "Nrhs", $Failed];
        Nlhs = Lookup[Eqs, "Nlhs", $Failed];
        bool = (EqEntryIn&&EqValueAuxiliaryEdges&&AltOptCond&&
            EqBalanceSplittingFlows&&AltFlows&&AltTransitionFlows&&
            IneqJs&&IneqJts&&IneqSwitchingByVertex)/.assoc;
        If[ bool,
            Print["All restrictions are ", styleg@"True"],
            Print["At least one of the restrictions is ", styler@bool];
            Print[style@"EqEntryIn: ", EqEntryIn/.assoc, "\n", EqEntryIn];
            Print[style@"EqValueAuxiliaryEdges: ", EqValueAuxiliaryEdges/.assoc, "\n", EqValueAuxiliaryEdges];
            Print[style@"AltOptCond: ", AltOptCond/.assoc, "\n", AltOptCond];
            Print[style@"EqBalanceSplittingFlows: ", EqBalanceSplittingFlows/.assoc, "\n", EqBalanceSplittingFlows];
            Print[style@"AltFlows: ", AltFlows/.assoc, "\n",AltFlows];
            Print[style@"AltTransitionFlows: ", AltTransitionFlows/.assoc, "\n", AltTransitionFlows];
            Print[style@"IneqJs: ", IneqJs/.assoc, "\n", IneqJs];
            Print[style@"IneqJts: ", IneqJts/.assoc, "\n", IneqJts];
            Print[style@"IneqSwitchingByVertex: ", IneqSwitchingByVertex/.assoc, "\n", IneqSwitchingByVertex]
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
    -1/m^2

A = 0.5;

(*V::usage = 
""
V = Function[{x, edge}, W[x, A]];

W::usage =
""
W = Function[{y, a}, a Sin[2 Pi (y + 1/4)]^2];*)

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
        jay = Eqs["SignedFlows"][List@@edge]/.sol;
        uT = u@@Reverse@(List@@edge);
        (*Print[uT];*)
        uT = uT/.sol;
        (*Print[uT];*)
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

plotM[Eqs_Association, string_String, pair_List] :=
    Module[ {jays, sol, edge},
        sol = Eqs[string];
        jays = Eqs["SignedFlows"][pair] /. sol;
        edge = UndirectedEdge@@pair;
        Plot[M[jays /. sol, x, edge], {x, 0, 1}, 
         PlotLabel -> edge,(*PlotRange\[Rule]{-0.1,2.4},*)
         GridLines -> Automatic]
    ];
   
plotMs[Eqs_, string_] :=
    plotM[Eqs, string, #] & /@ Eqs["EL"]
    
plotU[Eqs_Association, string_String, pair_List] :=
    Module[ {sol, edge},
        sol = Eqs[string];
        edge = UndirectedEdge@@pair;
        Plot[U[x, edge, Eqs, sol], {x, 0, 1}, PlotLabel -> edge(*,
         PlotRange\[Rule]{1-0.1,4.7}*), GridLines -> Automatic]
    ];

plotUs[Eqs_, string_] :=
    plotU[Eqs, string, #] & /@ Eqs["BEL"]