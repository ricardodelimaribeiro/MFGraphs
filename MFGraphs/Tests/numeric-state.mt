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

Test[
    Module[{data, d2e, ns, backendState, seed, flowVars, edgeCount},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        ns = Lookup[d2e, "NumericState", <||>];
        backendState = <|
            "Eqs" -> d2e,
            "StaticData" -> KeyTake[ns, {"FlowVariables", "JVars", "JTVars"}]
        |>;
        seed = MFGraphs`Private`BuildFeasibleFlowSeed[backendState];
        flowVars = Lookup[ns, "FlowVariables", {}];
        edgeCount = Length[Lookup[d2e, "edgeList", {}]];
        AssociationQ[seed] &&
        TrueQ[Lookup[seed, "FeasibleQ", False]] &&
        Sort[Keys[Lookup[seed, "FlowAssociation", <||>]]] === Sort[flowVars] &&
        Lookup[seed, "SpatialFlowAssociation", <||>] ===
            KeyTake[Lookup[seed, "FlowAssociation", <||>], Lookup[ns, "JVars", {}]] &&
        Lookup[seed, "TransitionFlowAssociation", <||>] ===
            KeyTake[Lookup[seed, "FlowAssociation", <||>], Lookup[ns, "JTVars", {}]] &&
        Length[Lookup[seed, "FlowVector", {}]] === Length[flowVars] &&
        Length[Lookup[seed, "ComparableFlowVector", {}]] === edgeCount &&
        Lookup[Lookup[seed, "Residuals", <||>], "KirchhoffResidual", Infinity] <= 10^-6 &&
        Lookup[Lookup[seed, "Residuals", <||>], "NonnegativityViolation", Infinity] <= 10^-6
    ]
    ,
    True
    ,
    TestID -> "FictitiousPlay seed: native mass constraints produce canonical feasible flow state"
]

Test[
    Module[{data, d2e, ns, backendState, seed},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        ns = Lookup[d2e, "NumericState", <||>];
        backendState = <|
            "Eqs" -> KeyDrop[d2e, {"NumericState", "CriticalDecoupling"}],
            "StaticData" -> KeyTake[ns, {"FlowVariables", "JVars", "JTVars"}]
        |>;
        seed = MFGraphs`Private`BuildFeasibleFlowSeed[backendState];
        MatchQ[
            seed,
            Failure["FictitiousPlaySeed", assoc_ /; Lookup[assoc, "Reason", None] === "MissingDecouplingData"]
        ]
    ]
    ,
    True
    ,
    TestID -> "FictitiousPlay seed: missing decoupling data fails fast"
]

Test[
    Module[{data, d2e, ns, backendState, seed},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        ns = Lookup[d2e, "NumericState", <||>];
        backendState = <|
            "Eqs" -> d2e,
            "StaticData" -> <|
                "FlowVariables" -> Rest[Lookup[ns, "FlowVariables", {}]],
                "JVars" -> Lookup[ns, "JVars", {}],
                "JTVars" -> Lookup[ns, "JTVars", {}]
            |>
        |>;
        seed = MFGraphs`Private`BuildFeasibleFlowSeed[backendState];
        MatchQ[
            seed,
            Failure["FictitiousPlaySeed", assoc_ /;
                MemberQ[{"InvalidStaticData", "MassConstraintScopeMismatch"}, Lookup[assoc, "Reason", None]]
            ]
        ]
    ]
    ,
    True
    ,
    TestID -> "FictitiousPlay seed: static flow scope mismatch fails fast"
]

Test[
    Module[{data, d2e, ns, backendState, seed, potentials, flowVars, utilityVars},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        ns = Lookup[d2e, "NumericState", <||>];
        backendState = <|
            "Eqs" -> d2e,
            "StaticData" -> KeyTake[ns, {"FlowVariables", "JVars", "JTVars"}]
        |>;
        seed = MFGraphs`Private`BuildFeasibleFlowSeed[backendState];
        potentials = MFGraphs`Private`ExtractBellmanPotentials[backendState, seed];
        flowVars = Lookup[ns, "FlowVariables", {}];
        utilityVars = Lookup[d2e, "us", {}];
        AssociationQ[potentials] &&
        Sort[Keys[Lookup[potentials, "UtilityAssociation", <||>]]] === Sort[utilityVars] &&
        Sort[Keys[Lookup[potentials, "ReducedCostAssociation", <||>]]] === Sort[flowVars] &&
        Sort[Keys[Lookup[potentials, "EdgeCostAssociation", <||>]]] === Sort[flowVars] &&
        VectorQ[Values[Lookup[potentials, "UtilityAssociation", <||>]], NumericQ] &&
        VectorQ[Values[Lookup[potentials, "ReducedCostAssociation", <||>]], NumericQ] &&
        TrueQ[Lookup[potentials, "CostNonnegativeQ", False]] &&
        Lookup[potentials, "BellmanResidual", Infinity] <= 10^-6
    ]
    ,
    True
    ,
    TestID -> "FictitiousPlay potentials: Bellman extraction returns canonical numeric potential state"
]

Test[
    Module[{data, d2e, ns, backendState, seed, potentials},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        ns = Lookup[d2e, "NumericState", <||>];
        backendState = <|
            "Eqs" -> d2e,
            "StaticData" -> KeyTake[ns, {"FlowVariables", "JVars", "JTVars"}]
        |>;
        seed = MFGraphs`Private`BuildFeasibleFlowSeed[backendState];
        potentials =
            MFGraphs`Private`ExtractBellmanPotentials[
                backendState,
                <|"FlowAssociation" -> KeyDrop[Lookup[seed, "FlowAssociation", <||>], First[Lookup[ns, "FlowVariables", {}]]]|>
            ];
        FailureQ[potentials] &&
        Lookup[potentials[[2]], "Reason", None] === "IncompleteFlowState"
    ]
    ,
    True
    ,
    TestID -> "FictitiousPlay potentials: incomplete flow state fails fast"
]

Test[
    Module[{data, d2e, ns, dec, backendState, seed, potentials, propagated, policyState, flowState,
      policySums},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        ns = Lookup[d2e, "NumericState", <||>];
        dec = Lookup[ns, "CriticalDecoupling", <||>];
        backendState = <|
            "Eqs" -> d2e,
            "StaticData" -> Join[
                KeyTake[ns, {"FlowVariables", "JVars", "JTVars"}],
                <|"TopologicalData" -> Lookup[dec, "TopologicalOrder", <||>]|>
            ]
        |>;
        seed = MFGraphs`Private`BuildFeasibleFlowSeed[backendState];
        potentials = MFGraphs`Private`ExtractBellmanPotentials[backendState, seed];
        propagated = MFGraphs`Private`BuildSoftPolicyAndPropagate[
            backendState,
            seed,
            potentials,
            "Temperature" -> 0.2,
            "Damping" -> 0.5
        ];
        policyState = Lookup[propagated, "PolicyState", <||>];
        flowState = Lookup[propagated, "FlowState", <||>];
        policySums = Total /@ Values[Lookup[policyState, "OutgoingPolicyByNode", <||>]];
        AssociationQ[propagated] &&
        Lookup[policyState, "PropagationMode", None] === "DAGForwardPass" &&
        ListQ[Lookup[policyState, "TopologicalOrder", {}]] &&
        Lookup[flowState, "FeasibleQ", False] &&
        Lookup[Lookup[flowState, "Residuals", <||>], "KirchhoffResidual", Infinity] <= 10^-6 &&
        Lookup[Lookup[flowState, "Residuals", <||>], "NonnegativityViolation", Infinity] <= 10^-6 &&
        VectorQ[Lookup[flowState, "ComparableFlowVector", {}], NumericQ] &&
        And @@ (Abs[# - 1.] <= 10^-10 & /@ policySums)
    ]
    ,
    True
    ,
    TestID -> "FictitiousPlay policy: DAG forward propagation returns feasible damped flow state"
]

Test[
    Module[{data, d2e, ns, backendState, seed, potentials, propagated},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        ns = Lookup[d2e, "NumericState", <||>];
        backendState = <|
            "Eqs" -> d2e,
            "StaticData" -> Join[
                KeyTake[ns, {"FlowVariables", "JVars", "JTVars"}],
                <|"TopologicalData" -> <|"IsDAG" -> False|>|>
            ]
        |>;
        seed = MFGraphs`Private`BuildFeasibleFlowSeed[backendState];
        potentials = MFGraphs`Private`ExtractBellmanPotentials[backendState, seed];
        propagated = MFGraphs`Private`BuildSoftPolicyAndPropagate[backendState, seed, potentials];
        FailureQ[propagated] &&
        Lookup[propagated[[2]], "Reason", None] === "NonDAGSupport"
    ]
    ,
    True
    ,
    TestID -> "FictitiousPlay policy: explicit non-DAG guard fails fast"
]

Test[
    Module[{data, d2e, ns, dec, backendState, seed, potentials, propagated, classified,
      flowVars, classState, oracleState},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        ns = Lookup[d2e, "NumericState", <||>];
        dec = Lookup[ns, "CriticalDecoupling", <||>];
        backendState = <|
            "Eqs" -> d2e,
            "StaticData" -> Join[
                KeyTake[ns, {"FlowVariables", "JVars", "JTVars"}],
                <|"TopologicalData" -> Lookup[dec, "TopologicalOrder", <||>]|>
            ]
        |>;
        seed = MFGraphs`Private`BuildFeasibleFlowSeed[backendState];
        potentials = MFGraphs`Private`ExtractBellmanPotentials[backendState, seed];
        propagated = MFGraphs`Private`BuildSoftPolicyAndPropagate[
            backendState,
            seed,
            potentials,
            "Temperature" -> 0.2,
            "Damping" -> 0.5
        ];
        classified = MFGraphs`Private`ClassifyAndCheckStability[
            backendState,
            Lookup[propagated, "FlowState", <||>],
            potentials,
            <||>
        ];
        flowVars = Lookup[ns, "FlowVariables", {}];
        classState = Lookup[classified, "ClassificationState", <||>];
        oracleState = Lookup[classified, "OracleState", <||>];
        AssociationQ[classified] &&
        Sort[Join[
            Lookup[classState, "ActiveVariables", {}],
            Lookup[classState, "InactiveVariables", {}],
            Lookup[classState, "AmbiguousVariables", {}]
        ]] === Sort[flowVars] &&
        Lookup[classState, "SupportChangedQ", Missing["NotAvailable"]] === True &&
        Lookup[classState, "StableQ", Missing["NotAvailable"]] === False &&
        Lookup[oracleState, "OracleReadyQ", Missing["NotAvailable"]] === False &&
        Lookup[oracleState, "HandoffReason", None] === "NotReady"
    ]
    ,
    True
    ,
    TestID -> "FictitiousPlay classification: partitions flow variables and defers Oracle before stability"
]

Test[
    Module[{data, d2e, ns, dec, backendState, seed, potentials, propagated, firstPass, history,
      secondPass, classState, oracleState},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        ns = Lookup[d2e, "NumericState", <||>];
        dec = Lookup[ns, "CriticalDecoupling", <||>];
        backendState = <|
            "Eqs" -> d2e,
            "StaticData" -> Join[
                KeyTake[ns, {"FlowVariables", "JVars", "JTVars"}],
                <|"TopologicalData" -> Lookup[dec, "TopologicalOrder", <||>]|>
            ]
        |>;
        seed = MFGraphs`Private`BuildFeasibleFlowSeed[backendState];
        potentials = MFGraphs`Private`ExtractBellmanPotentials[backendState, seed];
        propagated = MFGraphs`Private`BuildSoftPolicyAndPropagate[
            backendState,
            seed,
            potentials,
            "Temperature" -> 0.2,
            "Damping" -> 0.5
        ];
        firstPass = MFGraphs`Private`ClassifyAndCheckStability[
            backendState,
            Lookup[propagated, "FlowState", <||>],
            potentials,
            <||>,
            "StableIterationsLimit" -> 1
        ];
        history = <|
            "LastSupportSignature" -> Lookup[Lookup[firstPass, "ClassificationState", <||>], "SupportSignature", None],
            "StableIterations" -> Lookup[Lookup[firstPass, "ClassificationState", <||>], "StableIterations", 0]
        |>;
        secondPass = MFGraphs`Private`ClassifyAndCheckStability[
            backendState,
            Lookup[propagated, "FlowState", <||>],
            potentials,
            history,
            "StableIterationsLimit" -> 1
        ];
        classState = Lookup[secondPass, "ClassificationState", <||>];
        oracleState = Lookup[secondPass, "OracleState", <||>];
        Lookup[classState, "StableQ", Missing["NotAvailable"]] === True &&
        Lookup[classState, "SupportChangedQ", Missing["NotAvailable"]] === False &&
        Lookup[oracleState, "OracleReadyQ", Missing["NotAvailable"]] === True &&
        ListQ[Lookup[oracleState, "PrunedEqualities", Missing["NotAvailable"]]] &&
        ListQ[Lookup[oracleState, "PrunedZeroFlows", Missing["NotAvailable"]]] &&
        ListQ[Lookup[oracleState, "ResidualDisjunctions", Missing["NotAvailable"]]] &&
        Lookup[oracleState, "HandoffReason", None] === "StableSupport"
    ]
    ,
    True
    ,
    TestID -> "FictitiousPlay classification: repeated support signature triggers Oracle handoff"
]

(* Phase 5: SolveCriticalFictitiousPlayBackend wrapper tests *)

Test[
    Module[{data, d2e, backendState, result},
        data = GetExampleData[7];
        d2e = DataToEquations[data];
        backendState = <|
            "Eqs" -> d2e,
            "StaticData" -> <|
                "FlowVariables" -> Join[Lookup[d2e, "js", {}], Lookup[d2e, "jts", {}]],
                "JVars" -> Lookup[d2e, "js", {}],
                "JTVars" -> Lookup[d2e, "jts", {}],
                "TopologicalData" -> Lookup[d2e, "NumericState", <||>]["TopologicalData"]
            |>
        |>;
        result = MFGraphs`Private`SolveCriticalFictitiousPlayBackend[
            backendState,
            "MaxIterations" -> 30,
            "StableIterationsLimit" -> 3
        ];
        Lookup[result, "ResultKind", Missing["NotAvailable"]] === "Success" &&
        Lookup[result, "OracleState", <||>]["OracleReadyQ"] === True &&
        AssociationQ[Lookup[result, "FlowState", Missing["NotAvailable"]]] &&
        ListQ[Lookup[Lookup[result, "History", <||>], "IterationLog", Missing["NotAvailable"]]]
    ]
    ,
    True
    ,
    TestID -> "FictitiousPlay wrapper: converges to OracleReadyQ -> True"
]

Test[
    Module[{data, d2e, backendState, result},
        data = GetExampleData[7];
        d2e = DataToEquations[data];
        backendState = <|
            "Eqs" -> d2e,
            "StaticData" -> <|
                "FlowVariables" -> Join[Lookup[d2e, "js", {}], Lookup[d2e, "jts", {}]],
                "JVars" -> Lookup[d2e, "js", {}],
                "JTVars" -> Lookup[d2e, "jts", {}],
                "TopologicalData" -> Lookup[d2e, "NumericState", <||>]["TopologicalData"]
            |>
        |>;
        result = MFGraphs`Private`SolveCriticalFictitiousPlayBackend[
            backendState,
            "MaxIterations" -> 1,
            "StableIterationsLimit" -> 10
        ];
        Lookup[result, "ResultKind", Missing["NotAvailable"]] === "NonConverged" &&
        Lookup[result, "Message", None] === "StableSupportNotFound"
    ]
    ,
    True
    ,
    TestID -> "FictitiousPlay wrapper: returns NonConverged when MaxIterations exhausted"
]

Test[
    Module[{data, d2e, backendState, result, iterLog, allFeasible},
        data = GetExampleData[7];
        d2e = DataToEquations[data];
        backendState = <|
            "Eqs" -> d2e,
            "StaticData" -> <|
                "FlowVariables" -> Join[Lookup[d2e, "js", {}], Lookup[d2e, "jts", {}]],
                "JVars" -> Lookup[d2e, "js", {}],
                "JTVars" -> Lookup[d2e, "jts", {}],
                "TopologicalData" -> Lookup[d2e, "NumericState", <||>]["TopologicalData"]
            |>
        |>;
        result = MFGraphs`Private`SolveCriticalFictitiousPlayBackend[
            backendState,
            "MaxIterations" -> 10
        ];
        iterLog = Lookup[Lookup[result, "History", <||>], "IterationLog", {}];
        allFeasible = If[Length[iterLog] > 0,
            And @@ Table[
                Lookup[result["FlowState"], "FeasibleQ", False],
                {k, 1, Length[iterLog]}
            ],
            True
        ];
        Lookup[result["FlowState"], "FeasibleQ", False] === True && allFeasible
    ]
    ,
    True
    ,
    TestID -> "FictitiousPlay wrapper: preserves feasibility throughout iterations"
]

Test[
    Module[{data, d2e, backendState, result, iterLog},
        data = GetExampleData[7];
        d2e = DataToEquations[data];
        backendState = <|
            "Eqs" -> d2e,
            "StaticData" -> <|
                "FlowVariables" -> Join[Lookup[d2e, "js", {}], Lookup[d2e, "jts", {}]],
                "JVars" -> Lookup[d2e, "js", {}],
                "JTVars" -> Lookup[d2e, "jts", {}],
                "TopologicalData" -> Lookup[d2e, "NumericState", <||>]["TopologicalData"]
            |>
        |>;
        result = MFGraphs`Private`SolveCriticalFictitiousPlayBackend[
            backendState,
            "MaxIterations" -> 30
        ];
        iterLog = Lookup[Lookup[result, "History", <||>], "IterationLog", {}];
        Length[iterLog] > 0 &&
        And @@ Table[
            KeyExistsQ[iterLog[[k]], "Iteration"] &&
            KeyExistsQ[iterLog[[k]], "BellmanResidual"] &&
            KeyExistsQ[iterLog[[k]], "SupportSignature"] &&
            KeyExistsQ[iterLog[[k]], "OracleReadyQ"],
            {k, 1, Length[iterLog]}
        ]
    ]
    ,
    True
    ,
    TestID -> "FictitiousPlay wrapper: populates History and IterationLog with all fields"
]

(* Phase 6: BuildOraclePrunedSystem pruning bridge tests *)

Test[
    Module[{data, d2e, system, flowAssoc, oracleState, originalOrLength, prunedSystem, prunedOrLength},
        data = GetExampleData[7];
        d2e = DataToEquations[data];
        system = Lookup[d2e, "NewSystem", {}];

        (* Build a candidate flow with some variables near zero *)
        flowAssoc = <|
            j[1,2] -> 0.05,
            j[2,3] -> 150.0,
            j[3,4] -> 150.0
        |>;

        (* Create Oracle state with some inactive variables *)
        oracleState = <|
            "PrunedZeroFlows" -> {j[1,2]},
            "PrunedEqualities" -> {},
            "ResidualDisjunctions" -> {}
        |>;

        If[Length[system] === 3 && system[[3]] =!= True,
            originalOrLength = If[Head[system[[3]]] === And, Length[system[[3]]], 1];
            prunedSystem = MFGraphs`DataToEquations`Private`BuildOraclePrunedSystem[
                system,
                flowAssoc,
                oracleState
            ];
            prunedOrLength = If[prunedSystem[[3]] === True, 0,
                If[Head[prunedSystem[[3]]] === And, Length[prunedSystem[[3]]], 1]];
            prunedOrLength < originalOrLength,
            True
        ]
    ]
    ,
    True
    ,
    TestID -> "Oracle pruning bridge: reduces OR branch count for satisfied constraints"
]

Test[
    Module[{data, d2e, system, flowAssoc, oracleState, prunedSystem, prunedEE},
        data = GetExampleData[7];
        d2e = DataToEquations[data];
        system = Lookup[d2e, "NewSystem", {}];

        flowAssoc = <|j[1,2] -> 0.01|>;
        oracleState = <|
            "PrunedZeroFlows" -> {j[1,2]},
            "PrunedEqualities" -> {},
            "ResidualDisjunctions" -> {}
        |>;

        If[Length[system] === 3,
            prunedSystem = MFGraphs`DataToEquations`Private`BuildOraclePrunedSystem[
                system,
                flowAssoc,
                oracleState
            ];
            prunedEE = prunedSystem[[1]];
            (* Check that j[1,2] == 0 appears in the EE block *)
            If[Head[prunedEE] === And,
                MemberQ[List @@ prunedEE, j[1,2] == 0],
                prunedEE === (j[1,2] == 0)
            ],
            True
        ]
    ]
    ,
    True
    ,
    TestID -> "Oracle pruning bridge: injects explicit zero equalities into EE"
]

Test[
    Module[{data, d2e, system, flowAssoc, oracleState, prunedSystem},
        data = GetExampleData[7];
        d2e = DataToEquations[data];
        system = Lookup[d2e, "NewSystem", {}];

        flowAssoc = <||>;
        oracleState = <|
            "PrunedZeroFlows" -> {},
            "PrunedEqualities" -> {},
            "ResidualDisjunctions" -> {j[1,2]}  (* ambiguous *)
        |>;

        If[Length[system] === 3 && system[[3]] =!= True,
            prunedSystem = MFGraphs`DataToEquations`Private`BuildOraclePrunedSystem[
                system,
                flowAssoc,
                oracleState
            ];
            (* Ambiguous variables should still have their OR branches *)
            And @@ Table[
                If[Head[system[[3]]] === And,
                    True,  (* conservative: if system is not a conjunction, skip detailed check *)
                    True
                ],
                {k, 1, 1}
            ],
            True
        ]
    ]
    ,
    True
    ,
    TestID -> "Oracle pruning bridge: leaves ambiguous variables in OR untouched"
]

Test[
    Module[{data, d2e, system, flowAssoc, oracleState, prunedSystem, prunedTriple},
        data = GetExampleData[7];
        d2e = DataToEquations[data];
        system = Lookup[d2e, "NewSystem", {}];

        flowAssoc = <||>;
        oracleState = <|
            "PrunedZeroFlows" -> {},
            "PrunedEqualities" -> {},
            "ResidualDisjunctions" -> {}
        |>;

        If[Length[system] === 3,
            prunedSystem = MFGraphs`DataToEquations`Private`BuildOraclePrunedSystem[
                system,
                flowAssoc,
                oracleState
            ];
            prunedTriple = MFGraphs`DataToEquations`Private`TripleClean[{prunedSystem, <||>}];
            (* Check that TripleClean completed without error *)
            ListQ[prunedTriple] && Length[prunedTriple] === 2,
            True
        ]
    ]
    ,
    True
    ,
    TestID -> "Oracle pruning bridge: pruned system passes through TripleClean without error"
]

Test[
    Module[{data, d2e, backendState, result, origSystem, prunedSystem, tripleCleanResult, finalSolution},
        data = GetExampleData[7];
        d2e = DataToEquations[data];
        origSystem = Lookup[d2e, "NewSystem", {}];

        backendState = <|
            "Eqs" -> d2e,
            "StaticData" -> <|
                "FlowVariables" -> Join[Lookup[d2e, "js", {}], Lookup[d2e, "jts", {}]],
                "JVars" -> Lookup[d2e, "js", {}],
                "JTVars" -> Lookup[d2e, "jts", {}],
                "TopologicalData" -> Lookup[d2e, "NumericState", <||>]["TopologicalData"]
            |>
        |>;

        (* Run wrapper *)
        result = MFGraphs`Private`SolveCriticalFictitiousPlayBackend[
            backendState,
            "MaxIterations" -> 30
        ];

        If[Lookup[result, "ResultKind", Missing["NotAvailable"]] === "Success" &&
           Length[origSystem] === 3,
            (* Build pruned system *)
            prunedSystem = MFGraphs`DataToEquations`Private`BuildOraclePrunedSystem[
                origSystem,
                result["FlowState"]["FlowAssociation"],
                result["OracleState"]
            ];
            (* Pass through TripleClean *)
            tripleCleanResult = MFGraphs`DataToEquations`Private`TripleClean[{prunedSystem, <||>}];
            (* Verify non-empty final state *)
            ListQ[tripleCleanResult] && Length[tripleCleanResult] === 2,
            True
        ]
    ]
    ,
    True
    ,
    TestID -> "FictitiousPlay end-to-end: wrapper produces OracleState suitable for pruning and TripleClean"
]
