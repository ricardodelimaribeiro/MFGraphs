#!/usr/bin/env wolframscript
(* Usage: wolframscript -file Scripts/CompareSolverSolutions.wls [solver1] [solver2] *)
(* Compares residuals and rules from stored benchmark results using Boolean normalization. *)

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];

$Solvers = {"dnf", "optimizeddnf", "activeset", "reduce"};

Module[{args = Rest[$ScriptCommandLine], s1, s2, data1, data2, cases, commonCases},
  s1 = If[Length[args] >= 1, args[[1]], "dnf"];
  s2 = If[Length[args] >= 2, args[[2]], "optimizeddnf"];
  
  If[!MemberQ[$Solvers, s1] || !MemberQ[$Solvers, s2],
    Print["Error: Solvers must be one of: ", StringRiffle[$Solvers, ", "]];
    Exit[1]
  ];

  Print["Comparing Solver Solutions: ", s1, " vs ", s2];
  Print["=========================================="];

  data1 = Import[FileNameJoin[{$RepoRoot, "Results", "system_solver_" <> s1 <> "_solutions_latest.wl"}]];
  data2 = Import[FileNameJoin[{$RepoRoot, "Results", "system_solver_" <> s2 <> "_solutions_latest.wl"}]];

  If[!AssociationQ[data1] || !AssociationQ[data2],
    Print["Error: Could not load latest solution data for one or both solvers."];
    Exit[1]
  ];

  cases = Union[Keys[data1["Cases"]], Keys[data2["Cases"]]];
  commonCases = Intersection[Keys[data1["Cases"]], Keys[data2["Cases"]]];

  NormalizeResult[sol_] := Module[{rules, eqs},
    If[ListQ[sol],
      {Sort[sol], True},
      {Sort[Lookup[sol, "Rules", {}]], Lookup[sol, "Equations", True]}
    ]
  ];

  EquivalentQ[eq1_, eq2_] := Module[{b1, b2},
    If[eq1 === eq2, Return[True]];
    b1 = Quiet @ BooleanConvert[eq1, "DNF"];
    b2 = Quiet @ BooleanConvert[eq2, "DNF"];
    If[b1 === b2, Return[True]];
    (* Fallback to logical equivalence check if DNF forms differ syntactically *)
    Quiet @ Resolve[ForAll[{}, b1 == b2]] === True
  ];

  Print[StringForm["| Case | `` Rules | `` Residual | Equiv? |", s1, s2]];
  Print["|------|-----------|--------------|--------|"];

  Do[
    Module[{c1, c2, r1, e1, r2, e2, ruleMatch, eqMatch},
      c1 = data1["Cases"][name];
      c2 = data2["Cases"][name];
      
      If[MissingQ[c1] || MissingQ[c2] || c1["Status"] =!= "OK" || c2["Status"] =!= "OK",
        Print["| ", name, " | - | - | SKIP |"];
        Continue[]
      ];

      {r1, e1} = NormalizeResult[c1["Solution"]];
      {r2, e2} = NormalizeResult[c2["Solution"]];

      ruleMatch = (r1 === r2);
      eqMatch = EquivalentQ[e1, e2];

      Print["| ", name, " | ", 
        If[ruleMatch, "Match", "DIFF"], " | ",
        If[eqMatch, "Match", "DIFF"], " | ",
        If[ruleMatch && eqMatch, "YES", "NO"], " |"
      ];
      
      If[!(ruleMatch && eqMatch),
        If[!ruleMatch,
           Print["  [DIFF] Rules differ for ", name];
           (* Optional: Print diffs here if requested *)
        ];
        If[!eqMatch,
           Print["  [DIFF] Residuals differ for ", name];
           Print["    ", s1, ": ", Short[e1, 3]];
           Print["    ", s2, ": ", Short[e2, 3]];
        ];
      ];
    ],
    {name, commonCases}
  ];
];
