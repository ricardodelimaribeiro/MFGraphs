(* Exact verification against original generated constraint blocks. This module
   deliberately does not call buildSolverInputs, a package reducer, or an oracle. *)
BeginPackage["validationTools`", {"primitives`", "utilities`", "systemTools`"}];

exactSolutionReport::usage =
"exactSolutionReport[sys, sol] verifies a rule list or rules/residual family against \
all original generated constraints over the reals. Reports exact input provenance, \
nonemptiness, per-block soundness, parameter domain, and independently checks \
completeness by searching for an original solution outside the returned set. \
External family parameters are existentially projected for completeness; non-critical \
systems are reported as UnsupportedModel. \
Numerical inputs/results, unresolved conditions, infeasibility, and timeouts remain \
distinct. Options: \"TimeLimit\" (10 seconds for each of soundness and completeness), \
\"CheckCompleteness\" (True). No floating-point tolerance or rationalization is used.";

Begin["`Private`"];

$constraintNames = {"EqEntryIn", "EqBalanceSplittingFlows", "EqBalanceGatheringFlows",
    "EqGeneral", "IneqJs", "IneqJts", "IneqSwitchingByVertex", "IneqExitValues",
    "AltFlows", "AltTransitionFlows", "AltOptCond", "AltExitCond"};

originalConstraintBlocks[sys_] := AssociationMap[
    Function[name, With[{value = systemData[sys, name]},
        If[ListQ[value], And @@ value, value]]], $constraintNames];

(* Fresh scalar symbols make real quantification independent of WL's treatment
   of j[...] and u[...] as algebraic variables. *)
scalarVariables[expr_] := DeleteDuplicates @ Cases[expr,
    symbol_Symbol /; Context[symbol] =!= "System`", {0, Infinity}];

completedExpressionQ[expr_] := FreeQ[expr,
    _Reduce | _Resolve | _Solve | _FindInstance | _Failure | $Aborted | $TimedOut |
    $Failed | Indeterminate | _Missing | _ConditionalExpression | _Exists |
    _ForAll | _TerminatedEvaluation];

(* A successful backend must return an actual condition, not merely an expression
   absent from a list of known backend heads. *)
resolvedBooleanQ[True | False] := True;
resolvedBooleanQ[expr_And | expr_Or | expr_Not] := AllTrue[List @@ expr, resolvedBooleanQ];
resolvedBooleanQ[expr_] := completedExpressionQ[expr] && MatchQ[expr,
    _Equal | _Unequal | _Less | _LessEqual | _Greater | _GreaterEqual |
    _Inequality | _Element | _NotElement];

(* Complex numbers are atomic: searching only for _Real misses machine and
   arbitrary-precision Complex values. *)
inexactNumbers[expr_] := DeleteDuplicates @ Cases[expr,
    value_?NumberQ /; !ExactNumberQ[value], {0, Infinity}];

realReduce[expr_, vars_] := Quiet[Reduce[expr, vars, Reals]];
proofStatus[expr_] := Which[expr === False, "Proved", expr === $TimedOut,
    "Timeout", !resolvedBooleanQ[expr], "Unresolved", True, "Disproved"];

(* Independent reference path: Reduce the original linear equalities over Reals
   before the complement query. It shares neither the package's RowReduce
   substitution code nor its DNF engine. On a conjunction of exact affine
   equations Reduce's single branch is a complete parameterization. *)
independentComplement[original_, candidate_, vars_List] := Module[
    {equations, affine, rules, remaining, query, dnf, branches, bound},
    equations = Select[If[Head[original] === And, List @@ original, {original}], MatchQ[#, _Equal] &];
    affine = realReduce[And @@ equations, vars];
    If[affine === False, Return[False]];
    If[!resolvedBooleanQ[affine] || !FreeQ[affine, Or],
        Return[realReduce[original && Not[candidate], vars]]];
    rules = ToRules[affine];
    If[!ListQ[rules] || !AllTrue[rules, MatchQ[#, _Rule] &],
        Return[realReduce[original && Not[candidate], vars]]];
    query = (original && Not[candidate]) /. rules;
    remaining = Select[vars, !FreeQ[query, #] &];
    (* Independent complete Boolean expansion is economical on small systems.
       Reduce each conjunction; never discard a branch, timeout or unknown. *)
    bound = Times @@ Cases[query, expr_Or :> Length[expr], {0, Infinity}];
    If[bound <= 4096,
        dnf = BooleanConvert[query, "DNF"];
        branches = If[Head[dnf] === Or, List @@ dnf, {dnf}];
        Or @@ (realReduce[#, remaining] & /@ branches),
        realReduce[query, remaining]]
];

Options[exactSolutionReport] = {"TimeLimit" -> 10, "CheckCompleteness" -> True};

exactSolutionReport[sys_?mfgSystemQ, sol_, OptionsPattern[]] := Module[
    {limit = OptionValue["TimeLimit"], blocks, vars, scalar, toScalar, fromScalar,
     rules, residual, candidate, original, inputExact, resultExact, inexactInput,
     inexactResult, resolvedValues, normalizedRules, domain, free, domainReduced,
     nonempty = "NotChecked", blockChecks = <||>, soundness = "NotChecked",
     completeness = "NotChecked", complement = Missing["NotChecked"],
     kind = "Unresolved", checks, terminal, report, originalScalar,
     candidateScalar, checkTime, completenessTime = 0, scalarRules, allScalar,
     candidateEmpty = False, externalParameters, projectedCandidate, inputProvenance,
     families, realImage, nonemptyProof, dependencies},
    report[] := <|"ResultType" -> kind, "ExactInput" -> inputExact,
        "ExactOutput" -> resultExact, "InexactInputValues" -> inexactInput,
        "SourceProvenance" -> If[AssociationQ[inputProvenance], "Recorded", "LegacySystemMetadata"],
        "InexactOutputValues" -> inexactResult, "Nonempty" -> nonempty,
        "Soundness" -> soundness, "Completeness" -> completeness,
        "BlockChecks" -> blockChecks, "ParameterDomain" -> domainReduced,
        "FreeVariables" -> (free /. fromScalar),
        "SolutionExpression" -> candidate, "OriginalConstraints" -> original,
        "CompletenessCounterexampleSet" -> (complement /. fromScalar),
        "SoundnessSeconds" -> checkTime, "CompletenessSeconds" -> completenessTime|>;
    blocks = originalConstraintBlocks[sys];
    families = systemData[sys, #] & /@ {"Us", "Js", "Jts"};
    vars = If[AllTrue[families, ListQ], DeleteDuplicates[Join @@ families], {}];
    original = And @@ Values[blocks];
    inputProvenance = systemData[sys, "InputProvenance"];
    inexactInput = inexactNumbers[
        {original, systemData[sys, "Entries"], systemData[sys, "Exits"],
         systemData[sys, "SwitchingCosts"], inputProvenance}];
    inexactResult = inexactNumbers[sol];
    inputExact = inexactInput === {}; resultExact = inexactResult === {};
    fromScalar = {}; free = {}; checkTime = 0; domainReduced = Missing["NotChecked"];
    candidate = Missing["NotAvailable"];
    terminal = Which[sol === $TimedOut, "Timeout", sol === $Aborted, "Aborted",
        FailureQ[sol], "Failure", !FreeQ[blocks, _Missing] || !AllTrue[families, ListQ] ||
            !SubsetQ[vars, DeleteDuplicates @ Cases[original, _j | _u, Infinity]], "InvalidSystem",
        !criticalCongestionSystemQ[sys], "UnsupportedModel",
        !inputExact || !resultExact, "NumericalApproximation", True, None];
    If[terminal =!= None, kind = terminal; Return[report[]]];
    rules = Which[ListQ[sol], sol, AssociationQ[sol], Lookup[sol, "Rules", $Failed], True, $Failed];
    residual = If[AssociationQ[sol], Lookup[sol, "Residual", True], True];
    If[!ListQ[rules] || !AllTrue[rules, MatchQ[#, _Rule] && MemberQ[vars, First[#]] &] ||
        !DuplicateFreeQ[First /@ rules] || !completedExpressionQ[rules] ||
        !resolvedBooleanQ[residual], Return[report[]]];
    candidate = And[And @@ (Equal @@@ rules), residual];
    scalar = Table[Unique["exactVar$"], {Length[vars]}];
    toScalar = Thread[vars -> scalar]; fromScalar = Thread[scalar -> vars];
    originalScalar = original /. toScalar; candidateScalar = candidate /. toScalar;
    allScalar = DeleteDuplicates @ Join[scalar, scalarVariables[candidateScalar]];
    externalParameters = Complement[allScalar, scalar];
    If[candidate === False,
        kind = "InfeasibilityClaim"; nonempty = False; soundness = "Vacuous";
        If[TrueQ[OptionValue["CheckCompleteness"]],
            {completenessTime, complement} = AbsoluteTiming[TimeConstrained[
                realReduce[originalScalar, allScalar], limit, $TimedOut]];
            completeness = proofStatus[complement];
            If[completeness === "Proved", kind = "ExactInfeasible"]];
        Return[report[]]];
    scalarRules = rules /. toScalar;
    (* Reject dependencies that would grow or oscillate before attempting bounded
       substitution; self-identity rules impose no restriction and are harmless. *)
    dependencies = Flatten[Function[rule,
        DirectedEdge[First[rule], #] & /@ Select[scalar, !FreeQ[Last[rule], #] &]] /@
        Select[scalarRules, First[#] =!= Last[#] &]];
    If[!AcyclicGraphQ[Graph[scalar, dependencies]], Return[report[]]];
    resolvedValues = FixedPoint[(# /. scalarRules) &, scalar, Length[vars] + 1];
    (* Cyclic/nonterminating replacements cannot be used as a substitution proof. *)
    If[(resolvedValues /. scalarRules) =!= resolvedValues, Return[report[]]];
    normalizedRules = Thread[scalar -> resolvedValues];
    (* Substitution can cancel complex or undefined values in every model block.
       Keep the real-valued image obligation before using substitution as proof. *)
    realImage = And @@ (Element[#, Reals] & /@ resolvedValues);
    domain = (candidateScalar /. normalizedRules) && realImage;
    free = DeleteDuplicates @ Join[scalarVariables[resolvedValues], scalarVariables[domain]];
    {checkTime, checks} = AbsoluteTiming[TimeConstrained[
        domainReduced = realReduce[domain, free];
        If[resolvedBooleanQ[domainReduced],
            nonemptyProof = If[free === {}, domainReduced,
                With[{parameters = free, expression = domainReduced},
                    Quiet[Resolve[Exists[parameters, expression], Reals]]]];
            nonempty = Which[nonemptyProof === True, True, nonemptyProof === False, False,
                True, "Unresolved"];
            If[TrueQ[nonempty],
                blockChecks = Association @ KeyValueMap[
                    Function[{name, expr}, name -> proofStatus[realReduce[
                        domainReduced && Not[expr /. toScalar /. normalizedRules], free]]], blocks];
                soundness = If[AllTrue[Values[blockChecks], # === "Proved" &], "Proved",
                    If[MemberQ[Values[blockChecks], "Disproved"], "Disproved", "Unresolved"]],
                If[nonempty === False, candidateEmpty = True; soundness = "Vacuous",
                    soundness = "Unresolved"]],
            nonempty = "Unresolved"; soundness = "Unresolved"],
        limit, $TimedOut]];
    If[checks === $TimedOut, soundness = "Timeout"];
    domainReduced = domainReduced /. fromScalar;
    If[soundness === "Proved" && TrueQ[nonempty],
        kind = If[free === {}, "ExactDetermined", "ExactParametricFamily"]];
    If[soundness === "Disproved", kind = "InvalidSolution"];
    If[candidateEmpty, kind = "EmptyReturnedSet"];
    If[TrueQ[OptionValue["CheckCompleteness"]] && MemberQ[{"Proved", "Vacuous"}, soundness],
        (* This is an independent real quantifier-elimination query on the
           ORIGINAL blocks. No solver preprocessing, pruning, or equality cache. *)
        {completenessTime, complement} = AbsoluteTiming[TimeConstrained[
            If[candidateEmpty,
                (* The returned real set is already proved empty.  Its
                   complement within the original problem is simply the
                   original solution set; do not feed nonreal candidate
                   constants back into Reduce[..., Reals]. *)
                realReduce[originalScalar, scalar],
                projectedCandidate = If[externalParameters === {}, candidateScalar,
                    With[{parameters = externalParameters, expression = candidateScalar},
                        realReduce[Exists[parameters, expression], scalar]]];
                If[resolvedBooleanQ[projectedCandidate],
                    independentComplement[originalScalar, projectedCandidate, scalar],
                    projectedCandidate]
            ],
            limit, $TimedOut]];
        completeness = proofStatus[complement]];
    report[]
];

End[];
EndPackage[];
