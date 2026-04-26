(* Wolfram Language Package *)
(*
   MFGraphs: A Wolfram Language package for Mean Field Games on Networks.

   Thin loader — mirrors the structure of Maydan.m in the numerics package.

   PrependTo adds the MFGraphs directory to $Path so that WL can find each
   sub-package by context name. BeginPackage declares all sub-packages as
   dependencies; EndPackage appends all their contexts to $ContextPath, making
   every public symbol accessible unqualified after Needs["MFGraphs`"].

   Package loading DAG:
     primitives  (shared primitives: j, u, z, alpha, Cost, flags)
       ↓
     scenarioTools  (typed scenario kernel)
       ↓
     examples       (named scenario factories)
     unknownsTools  (symbolic unknown families)
       ↓
     systemTools    (structural equation system)
       ↓
     solver         (ReduceSystem solver)
     graphics       (visualization helpers)
*)

PrependTo[$Path, DirectoryName[$InputFileName]];

BeginPackage["MFGraphs`",
  { "primitives`",
    "scenarioTools`",
    "examples`",
    "unknownsTools`",
    "systemTools`",
    "solver`",
    "graphics`"
  }
];

EndPackage[]
