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
    Module[{data, d2e, ns, critical, assoc, vec, signedFast, residualFast, comp, diffSigned,
      entryPairs, exitPairs, boundaryInExpected, boundaryOutExpected},
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
        entryPairs = List @@@ Lookup[d2e, "entryEdges", {}];
        exitPairs = List @@@ Lookup[d2e, "exitEdges", {}];
        boundaryInExpected = N @ Total[Lookup[assoc, MFGraphs`Private`j @@@ entryPairs, 0]];
        boundaryOutExpected = N @ Total[Lookup[assoc, MFGraphs`Private`j @@@ exitPairs, 0]];
        VectorQ[signedFast, NumericQ] &&
        NumericQ[residualFast] &&
        NumericQ[Lookup[comp, "KirchhoffResidual", Missing["NotAvailable"]]] &&
        diffSigned <= 10^-10 &&
        Abs[residualFast - Lookup[comp, "KirchhoffResidual", 0.]] <= 10^-10 &&
        NumericQ[Lookup[comp, "BoundaryMassResidual", Missing["NotAvailable"]]] &&
        Abs[Lookup[comp, "BoundaryMassIn", Infinity] - boundaryInExpected] <= 10^-10 &&
        Abs[Lookup[comp, "BoundaryMassOut", Infinity] - boundaryOutExpected] <= 10^-10 &&
        Lookup[comp, "BoundaryMassResidual", Infinity] <= 10^-10
    ]
    ,
    True
    ,
    TestID -> "NumericState helpers: fast signed-flow, residual, and boundary-mass parity"
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

Test[
    Module[{data, d2e, pre, reduction, result, classResidual},
        data = GetExampleData[14] /. {I1 -> 10, U1 -> 0, S1 -> 0, S2 -> 0};
        d2e = DataToEquations[data];
        pre = MFGPreprocessing[d2e];
        reduction = Lookup[pre, "UtilityReduction", <||>];
        result = Quiet[CriticalCongestionSolver[pre]];
        classResidual = Lookup[result, "UtilityClassResidual", Missing["NotAvailable"]];
        AssociationQ[reduction] &&
        TrueQ[Lookup[reduction, "Enabled", False]] &&
        IntegerQ[Lookup[reduction, "ClassCount", Missing["NotAvailable"]]] &&
        Lookup[reduction, "ClassCount", 0] > 0 &&
        AssociationQ[Lookup[reduction, "RepresentativeRules", <||>]] &&
        NumericQ[classResidual] &&
        classResidual <= 10^-10
    ]
    ,
    True
    ,
    TestID -> "MFGPreprocessing utility reduction: zero-switching classes and residual telemetry"
]

Test[
    Module[{data, d2e, ns, decoupling},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        ns = Lookup[d2e, "NumericState", <||>];
        decoupling = Lookup[ns, "CriticalDecoupling", Missing["NotAvailable"]];
        AssociationQ[decoupling] &&
        And @@ (KeyExistsQ[decoupling, #] & /@ {"JVars", "UVars", "MassConstraints", "TopologicalOrder"})
    ]
    ,
    True
    ,
    TestID -> "NumericState decoupling: schema is attached"
]

Test[
    Module[{data, d2e, decoupling, massConstraints, aeq, beq, vars},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        decoupling = Lookup[d2e["NumericState"], "CriticalDecoupling", <||>];
        massConstraints = Lookup[decoupling, "MassConstraints", <||>];
        aeq = Lookup[massConstraints, "Aeq", Missing["NotAvailable"]];
        beq = Lookup[massConstraints, "beq", Missing["NotAvailable"]];
        vars = Lookup[massConstraints, "Variables", {}];
        Head[aeq] === SparseArray &&
        Developer`PackedArrayQ[beq] &&
        VectorQ[beq, NumericQ] &&
        Length[Dimensions[aeq]] === 2 &&
        Length[vars] === Dimensions[aeq][[2]] &&
        FreeQ[vars, _u]
    ]
    ,
    True
    ,
    TestID -> "NumericState decoupling: mass constraints use flow variables only"
]

Test[
    Module[{data, d2e1, d2e2, dec1, dec2, topo1, topo2, jVars, jVarOrder},
        data = GetExampleData[9] /. {I1 -> 100, I2 -> 40, U1 -> 0};
        d2e1 = DataToEquations[data];
        d2e2 = DataToEquations[data];
        dec1 = Lookup[d2e1["NumericState"], "CriticalDecoupling", <||>];
        dec2 = Lookup[d2e2["NumericState"], "CriticalDecoupling", <||>];
        topo1 = Lookup[dec1, "TopologicalOrder", <||>];
        topo2 = Lookup[dec2, "TopologicalOrder", <||>];
        jVars = Lookup[dec1, "JVars", {}];
        jVarOrder = Lookup[topo1, "JVarOrder", {}];
        TrueQ[Lookup[topo1, "IsDAG", False]] &&
        ListQ[Lookup[topo1, "NodeOrder", {}]] &&
        Length[Lookup[topo1, "NodeOrder", {}]] > 0 &&
        topo1 === topo2 &&
        And @@ (MemberQ[jVars, #] & /@ jVarOrder)
    ]
    ,
    True
    ,
    TestID -> "NumericState decoupling: DAG topological metadata is deterministic"
]

Test[
    Module[{data, d2e, dec, topo, cycleWitness},
        data = GetExampleData[27] /. {I1 -> 100, U1 -> 0};
        d2e = DataToEquations[data];
        dec = Lookup[d2e["NumericState"], "CriticalDecoupling", <||>];
        topo = Lookup[dec, "TopologicalOrder", <||>];
        cycleWitness = Lookup[topo, "CycleWitness", Missing["NotAvailable"]];
        TrueQ[Lookup[topo, "IsDAG", True] === False] &&
        Lookup[topo, "NodeOrder", {1}] === {} &&
        cycleWitness =!= Missing["NotApplicable"] &&
        cycleWitness =!= Missing["NotAvailable"] &&
        cycleWitness =!= {}
    ]
    ,
    True
    ,
    TestID -> "NumericState decoupling: non-DAG case exposes cycle witness"
]
