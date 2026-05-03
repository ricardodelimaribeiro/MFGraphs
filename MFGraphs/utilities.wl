(* Wolfram Language package *)
(* utilities.wl — shared logic for MFGraphs typed objects, solvers, and rules. *)

BeginPackage["utilities`", {"primitives`"}]

(* --- Public API declarations --- *)

mfgTypedQ::usage = "mfgTypedQ[obj, head] returns True if obj is of the form head[assoc_Association].";

mfgData::usage =
"mfgData[obj] returns the underlying Association of a typed object. \
mfgData[obj, key] returns the value for key, or Missing[\"KeyAbsent\", key].";

mergeRules::usage =
"mergeRules[oldRules, newRules] joins rule lists while keeping the latest rule for each left-hand side.";

normalizeRules::usage =
"normalizeRules[rules] rewrites right-hand sides through the full rule set.";

extractRules::usage =
"extractRules[sol] extracts replacement rules from a solution list or solution association.";

roundValues::usage = "roundValues[x] rounds numerical values in x to a standard precision (10^-10).";

exactBoundaryValue::usage =
"exactBoundaryValue[val] converts a numeric real boundary value to an exact rational or returns a Failure.";

exactBoundaryValues::usage =
"exactBoundaryValues[vals] applies exactBoundaryValue to a list and returns the first Failure if conversion fails.";

criticalCongestionSystemQ::usage =
"criticalCongestionSystemQ[sys] returns True if the system is a critical congestion system (Alpha == 1 everywhere).";

withCriticalCongestionSolver::usage =
"withCriticalCongestionSolver[sys, solverName, bodyFunc] handles validation, input building, \
and rule attachment for critical-congestion structural solvers.";

Begin["`Private`"]

(* --- Typed Object Infrastructure --- *)

mfgTypedQ[head_[_Association], head_] := True;
mfgTypedQ[_, _] := False;

mfgData[head_[assoc_Association]] := assoc;
mfgData[head_[assoc_Association], key_] := Lookup[assoc, key, Missing["KeyAbsent", key]];

(* --- Rule Management --- *)

mergeRules[oldRules_List, newRules_List] :=
    Reverse @ DeleteDuplicatesBy[Reverse @ Join[oldRules, newRules], First];

normalizeRules[rules_List] :=
    Map[Rule[First[#], ReplaceRepeated[Last[#], rules]] &, rules];

extractRules[sol_List] := Select[sol, MatchQ[#, _Rule | _RuleDelayed] &];
extractRules[sol_Association] := extractRules[Lookup[sol, "Rules", {}]];
extractRules[_] := {};

(* --- Numerics --- *)

roundValues[x_?NumberQ] := Round[x, 10^-10]
roundValues[Rule[a_, b_]] := Rule[a, roundValues[b]]
roundValues[x_List] := roundValues /@ x
roundValues[x_Association] := roundValues /@ x
roundValues[other_] := other

exactBoundaryValue[val_] /; ExactNumberQ[val] := val;
exactBoundaryValue[val_?NumericQ] :=
    Module[{nval = N[val]},
        If[TrueQ[Im[nval] == 0],
            Rationalize[val, 0],
            Failure["exactBoundaryValue", <|
                "Message" -> "Value must be a real numeric value.",
                "Value" -> val
            |>]
        ]
    ];
exactBoundaryValue[val_] :=
    Failure["exactBoundaryValue", <|
        "Message" -> "Value must be numeric.",
        "Value" -> val
    |>];

exactBoundaryValues[vals_List] :=
    Module[{exactVals = exactBoundaryValue /@ vals, firstFailure},
        firstFailure = FirstCase[exactVals, _Failure, Missing["NotFound"]];
        If[FailureQ[firstFailure], firstFailure, exactVals]
    ];

(* --- Solver Orchestration --- *)

(* Shared check for solver eligibility. *)
criticalCongestionSystemQ[sys_] :=
    Module[{halfPairs, hamiltonian, alphaDefault, edgeAlpha, alphaAtEdge},
        halfPairs   = mfgData[sys, "HalfPairs"];
        hamiltonian = mfgData[sys, "Hamiltonian"];
        If[MissingQ[hamiltonian] || !AssociationQ[hamiltonian], hamiltonian = <||>];
        alphaDefault = Lookup[hamiltonian, "Alpha", 1];
        edgeAlpha    = Lookup[hamiltonian, "EdgeAlpha", <||>];
        If[!ListQ[halfPairs] || !AssociationQ[edgeAlpha], Return[False, Module]];
        alphaAtEdge[edge_List] :=
            Lookup[
                edgeAlpha,
                Key[edge],
                Lookup[edgeAlpha, Key[Reverse[edge]], alphaDefault]
            ];
        alphaDefault === 1 && AllTrue[halfPairs, alphaAtEdge[#] === 1 &]
    ];

withCriticalCongestionSolver[sys_, solverName_String, bodyFunc_, buildInputsFunc_, attachRulesFunc_] :=
    Module[{inputs, constraints, allVars, rulesAcc, result},
        (* 1. Validation *)
        If[!criticalCongestionSystemQ[sys],
            With[{sym = Symbol[solverName]},
                Message[MessageName[sym, "noncritical"]]
            ];
            Return[Failure[solverName, <|"Message" -> solverName <> " supports only critical congestion systems with Alpha == 1 on every edge."|>], Module]
        ];

        (* 2. Preprocessing *)
        inputs = buildInputsFunc[sys];
        If[!ListQ[inputs] || Length[inputs] < 3,
            Return[Failure[solverName, <|"Message" -> "Failed to build solver inputs."|>], Module]
        ];
        {constraints, allVars, rulesAcc} = inputs;

        (* 3. Core execution *)
        result = bodyFunc[constraints, allVars];

        (* 4. Post-processing *)
        attachRulesFunc[result, rulesAcc]
    ];

End[]

EndPackage[]
