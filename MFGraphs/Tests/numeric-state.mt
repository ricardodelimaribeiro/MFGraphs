(* Wolfram Language Test file *)
(* Numeric state compiler + adapter contract tests (internal, no solver routing changes). *)

Test[
    Module[{data, d2e, ns},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        ns = Lookup[d2e, "NumericState", Missing["NotAvailable"]];
        AssociationQ[ns] &&
        Lookup[ns, "Version", None] === 1 &&
        Length[Lookup[ns, "FlowVariables", {}]] ===
            Length[Lookup[d2e, "js", {}]] + Length[Lookup[d2e, "jts", {}]] &&
        AssociationQ[Lookup[ns, "EdgeIndex", <||>]] &&
        AssociationQ[Lookup[ns, "TransitionIndex", <||>]] &&
        Head[Lookup[ns, "KirchhoffKM", None]] === SparseArray &&
        Head[Lookup[ns, "SignedEdgeMatrix", None]] === SparseArray
    ]
    ,
    True
    ,
    TestID -> "NumericState: schema is attached by DataToEquations"
]

Test[
    Module[{data, ns1, ns2},
        data = GetExampleData[12] /. {I1 -> 2, U1 -> 0};
        ns1 = DataToEquations[data]["NumericState"];
        ns2 = DataToEquations[data]["NumericState"];
        Lookup[ns1, {"EdgePairs", "TransitionTriples", "JVars", "JTVars", "UVars", "KirchhoffVariables"}] ===
            Lookup[ns2, {"EdgePairs", "TransitionTriples", "JVars", "JTVars", "UVars", "KirchhoffVariables"}] &&
        Normal[Lookup[ns1, "KirchhoffKM", SparseArray[{}]]] ===
            Normal[Lookup[ns2, "KirchhoffKM", SparseArray[{}]]] &&
        Normal[Lookup[ns1, "SignedEdgeMatrix", SparseArray[{}]]] ===
            Normal[Lookup[ns2, "SignedEdgeMatrix", SparseArray[{}]]] &&
        Lookup[ns1, "DistanceToExit", <||>] === Lookup[ns2, "DistanceToExit", <||>]
    ]
    ,
    True
    ,
    TestID -> "NumericState: deterministic compilation for identical input"
]

Test[
    Module[{data, d2e, ns, critical, assoc, vec},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        ns = d2e["NumericState"];
        critical = Quiet[CriticalCongestionSolver[d2e]];
        assoc = critical["AssoCritical"];
        vec = MFGraphs`Private`EncodeFlowAssociation[ns, assoc];
        Developer`PackedArrayQ[vec] &&
        VectorQ[vec, NumericQ] &&
        Length[vec] === Length[Lookup[ns, "FlowVariables", {}]]
    ]
    ,
    True
    ,
    TestID -> "NumericState adapters: encode returns packed numeric vector"
]

Test[
    Module[{data, d2e, ns, critical, assoc, vec, decoded, vars, maxDiff},
        data = GetExampleData[12] /. {I1 -> 2, U1 -> 0};
        d2e = DataToEquations[data];
        ns = d2e["NumericState"];
        critical = Quiet[CriticalCongestionSolver[d2e]];
        assoc = critical["AssoCritical"];
        vec = MFGraphs`Private`EncodeFlowAssociation[ns, assoc];
        decoded = MFGraphs`Private`DecodeFlowVector[ns, vec];
        vars = Lookup[ns, "FlowVariables", {}];
        maxDiff =
            Quiet @ Check[
                Max @ Abs @ N[Lookup[decoded, vars] - Lookup[assoc, vars]],
                Infinity
            ];
        AssociationQ[decoded] &&
        NumericQ[maxDiff] &&
        maxDiff <= 10^-10
    ]
    ,
    True
    ,
    TestID -> "NumericState adapters: encode/decode round-trip identity"
]

Test[
    Module[{data, d2e, ns, critical, assoc, vec, signedFast, residualFast, comp, diffSigned},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        ns = d2e["NumericState"];
        critical = Quiet[CriticalCongestionSolver[d2e]];
        assoc = critical["AssoCritical"];
        vec = MFGraphs`Private`EncodeFlowAssociation[ns, assoc];
        signedFast = MFGraphs`Private`ComputeSignedEdgeFlowsFast[ns, vec];
        residualFast = MFGraphs`Private`ComputeKirchhoffResidualFast[ns, vec];
        comp = MFGraphs`Private`BuildSolverComparisonData[d2e, assoc];
        diffSigned =
            Quiet @ Check[
                Norm[N[signedFast - Lookup[comp, "ComparableFlowVector", {}]], Infinity],
                Infinity
            ];
        VectorQ[signedFast, NumericQ] &&
        NumericQ[residualFast] &&
        NumericQ[Lookup[comp, "KirchhoffResidual", Missing["NotAvailable"]]] &&
        diffSigned <= 10^-10 &&
        Abs[residualFast - Lookup[comp, "KirchhoffResidual", 0.]] <= 10^-10
    ]
    ,
    True
    ,
    TestID -> "NumericState helpers: fast signed-flow and residual parity"
]

Test[
    Module[{decoded},
        decoded = MFGraphs`Private`DecodeFlowVector[
            <|"FlowVariables" -> {j[1, 2], j[2, 3]}|>,
            {1.0}
        ];
        decoded === <||>
    ]
    ,
    True
    ,
    TestID -> "NumericState adapters: decode guards mismatched vector length"
]
