#!/usr/bin/env wolframscript
(* ProfileMFGSS.wls -- Profile MFGSystemSolver internals step by step *)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

SetDirectory[DirectoryName[$InputFileName] // ParentDirectory];

Needs["MFGraphs`"];

Get["Scripts/BenchmarkHelpers.wls"];

$MFGraphsVerbose = False;

FmtTime[t_?NumericQ] := ToString[CForm[N[t]]];
FmtTime[x_] := ToString[x];

(* RoundValues is in MFGraphs`Private` — redefine locally *)
localRoundValues[x_?NumberQ] := Round[x, 10^-10];
localRoundValues[Rule[a_, b_]] := Rule[a, localRoundValues[b]];
localRoundValues[x_List] := localRoundValues /@ x;
localRoundValues[x_Association] := localRoundValues /@ x;
localRoundValues[x_] := x;

(* ============================================================ *)
(* Step-by-step MFGSystemSolver profiling                      *)
(* Manually replicates MFGSystemSolver pipeline with timing     *)
(* ============================================================ *)

ProfileMFGSS[caseKey_] := Module[
  {data, Eqs, PreEqs, js, jts, us, NewSystem, InitRules, costpluscurrents,
   approxJs, Ncpc, ineqs, ineqsByTransition, sorted,
   time, timings, dnfResult, reducedResult, System, vars,
   jjtsR, usR, pickOne},

  timings = <||>;

  (* Setup *)
  data = GetExampleData[caseKey] /. GetParams[caseKey];
  Eqs = DataToEquations[data];
  If[Eqs === $Failed, Print["  DataToEquations failed. Skipping."]; Return[<||>, Module]];

  ClearSolveCache[];
  {time, PreEqs} = AbsoluteTiming[MFGPreprocessing[Eqs]];
  timings["MFGPreprocessing"] = time;

  us = Lookup[PreEqs, "us"];
  js = Lookup[PreEqs, "js"];
  jts = Lookup[PreEqs, "jts"];
  InitRules = Lookup[PreEqs, "InitRules"];
  NewSystem = Lookup[PreEqs, "NewSystem"];
  costpluscurrents = Lookup[PreEqs, "costpluscurrents"];

  (* Simulate critical congestion: all flows = 0 *)
  approxJs = AssociationThread[js, 0 js];

  ClearSolveCache[];

  (* Step 1: Cost substitution *)
  {time, Ncpc} = AbsoluteTiming[localRoundValues @ (Expand /@ (costpluscurrents /. approxJs))];
  timings["1_CostSubstitution"] = time;

  {time, InitRules} = AbsoluteTiming[Expand /@ (InitRules /. Ncpc)];
  timings["2_ExpandInitRules"] = time;

  {time, NewSystem} = AbsoluteTiming[NewSystem /. Ncpc];
  timings["3_SubstituteNewSystem"] = time;

  (* Step 2: Group inequalities by transition flow *)
  If[NewSystem[[2]] === True,
    timings["4_GroupIneqByTransition"] = 0.;
    timings["5_SimplifyIneqByTransition"] = 0.;
    ineqsByTransition = ConstantArray[True, Length[jts]];
    ,
    ineqs = NewSystem[[2]];

    (* Baseline: Select approach *)
    {time, ineqsByTransition} = AbsoluteTiming[
      With[{iq = ineqs},
        MFGParallelMap[Function[jt, Select[iq, !FreeQ[jt][#]&]], jts]
      ]
    ];
    timings["4_GroupIneqByTransition_Select"] = time;

    (* Alternative: Index-based approach *)
    Module[{ineqList, index, timeIdx, resultIdx},
      ineqList = If[Head[ineqs] === And, List @@ ineqs, {ineqs}];
      {timeIdx, resultIdx} = AbsoluteTiming[
        index = Association[# -> {} & /@ jts];
        Do[
          Do[
            If[!FreeQ[ineq, jt], index[jt] = Append[index[jt], ineq]],
            {jt, jts}
          ],
          {ineq, ineqList}
        ];
        If[Length[#] == 0, True, And @@ #] & /@ Values[index]
      ];
      timings["4_GroupIneqByTransition_Index"] = timeIdx;
    ];

    (* Step 3: Simplify inequalities by transition flow *)
    {time, ineqsByTransition} = AbsoluteTiming[
      DeduplicateByComplexity[MFGParallelMap[Simplify, ineqsByTransition]]
    ];
    timings["5_SimplifyIneqByTransition"] = time;
  ];

  NewSystem[[2]] = (And @@ ineqsByTransition);
  NewSystem = SystemToTriple[And @@ NewSystem];

  (* Step 4: TripleClean *)
  {time, {NewSystem, InitRules}} = AbsoluteTiming[TripleClean[{NewSystem, InitRules}]];
  timings["6_TripleClean"] = time;

  (* Step 5: DNFSolveStep breakdown *)
  Module[{ns, nr, nsInner, nrInner, sortedInner, timeDNF, timeReduce, timeBC, timeTC2},
    {ns, nr} = TripleClean[{NewSystem, InitRules}];

    If[ns[[3]] === True,
      timings["7_DNFSolveStep_DNFReduce"] = 0.;
      timings["8_DNFSolveStep_ReduceDisjuncts"] = 0.;
      timings["9_DNFSolveStep_BooleanConvert"] = 0.;
      timings["10_DNFSolveStep_FinalTripleClean"] = 0.;
      ,
      sortedInner = DeduplicateByComplexity[ns[[3]]];

      Print["  DNF input: ", Length[sortedInner], " disjunctions, system size: ",
            If[Head[And @@ Most[ns]] === And, Length[And @@ Most[ns]], 1]];

      {timeDNF, nsInner} = AbsoluteTiming[DNFReduce[And @@ Most[ns], sortedInner]];
      timings["7_DNFSolveStep_DNFReduce"] = timeDNF;

      Print["  DNFReduce output: ",
            If[Head[nsInner] === Or, Length[nsInner], 1], " disjunct(s), LeafCount: ", LeafCount[nsInner]];

      {timeReduce, nsInner} = AbsoluteTiming[ReduceDisjuncts[nsInner]];
      timings["8_DNFSolveStep_ReduceDisjuncts"] = timeReduce;

      Print["  ReduceDisjuncts output: ",
            If[Head[nsInner] === Or, Length[nsInner], 1], " disjunct(s), LeafCount: ", LeafCount[nsInner]];

      {timeBC, nsInner} = AbsoluteTiming[BooleanConvert @ nsInner];
      timings["9_DNFSolveStep_BooleanConvert"] = timeBC;

      nsInner = SystemToTriple[nsInner];
      {timeTC2, {nsInner, nr}} = AbsoluteTiming[TripleClean[{nsInner, nr}]];
      timings["10_DNFSolveStep_FinalTripleClean"] = timeTC2;

      NewSystem = nsInner;
      InitRules = nr;
    ];
  ];

  (* Step 6: Final system simplification + FindInstance *)
  {time, System} = AbsoluteTiming[Simplify[And @@ NewSystem]];
  timings["11_FinalSimplify"] = time;

  {time, {System, InitRules}} = AbsoluteTiming[
    Module[{s, ir},
      {s, ir} = TripleClean[{SystemToTriple[System], InitRules}];
      {And @@ s, ir}
    ]
  ];
  timings["12_FinalTripleClean"] = time;

  If[System =!= True && System =!= False,
    vars = Select[Variables[System], MatchQ[#, j[_, _, _] | j[_, _] | u[_, _] | u[_, _, _]]&];
    vars = Complement[vars, Keys[InitRules]];
    jjtsR = Select[vars, MatchQ[#, j[_, _, _] | j[_, _]]&];

    If[Length[vars] > 0,
      {time, pickOne} = AbsoluteTiming[
        Module[{p},
          p = FindInstance[System && And @@ ((# > 0)& /@ jjtsR), vars, Reals];
          If[p === {},
            p = FindInstance[System && And @@ ((# >= 0)& /@ jjtsR), vars, Reals]
          ];
          p
        ]
      ];
      timings["13_FindInstance"] = time;
      ,
      timings["13_FindInstance"] = 0.;
    ];
    ,
    timings["13_FindInstance"] = 0.;
  ];

  (* Total *)
  timings["TOTAL"] = Total[Values[timings]];

  timings
];

(* ============================================================ *)
(* MAIN                                                         *)
(* ============================================================ *)

Print["===== MFGSystemSolver Detailed Profiler ====="];
Print[];

$TestCases = {8, 19, 14};

Do[
  Print["==== Case: ", caseKey, " ===="];
  Module[{timings},
    timings = ProfileMFGSS[caseKey];
    If[timings =!= <||>,
      Do[
        Print["  ", StringPadRight[key, 45], FmtTime[timings[key]], " s"],
        {key, Keys[timings]}
      ];
    ];
    Print[];
  ],
  {caseKey, $TestCases}
];

Print["Done."];

(* :!CodeAnalysis::EndBlock:: *)
