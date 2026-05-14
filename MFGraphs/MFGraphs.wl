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
     solversTools   (ReduceSystem solver)
     orchestrationTools (High-level orchestration)
     graphicsTools  (visualization helpers)

   Optional subpackages (NOT loaded here; load explicitly after Needs["MFGraphs`"]):
     Tawaf`         (unrolled circumambulation scenario builder)
     numericOracle` (LP-relaxation oracle; introduces FindInstance floats)
*)

PrependTo[$Path, DirectoryName[$InputFileName]];

BeginPackage["MFGraphs`",
  { "primitives`",
    "utilities`",
    "scenarioTools`",
    "examples`",
    "unknownsTools`",
    "systemTools`",
    "solversTools`",
    "orchestrationTools`",
    "graphicsTools`"
  }
];

EndPackage[]
