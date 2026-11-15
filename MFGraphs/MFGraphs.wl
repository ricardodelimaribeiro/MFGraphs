(* Wolfram Language Package *)
(*
   MFGraphs: A Wolfram Language package for Mean Field Games on Networks.

   Thin loader — mirrors the structure of Maydan.m in the numerics package.

   The first expression puts the MFGraphs directory at the front of $Path so
   that WL can find each sub-package by context name. It is load-bearing:
   PacletInfo.m declares only MFGraphs`, so the sub-contexts are found only
   through this $Path entry. Re-loading moves the entry back to the front
   instead of adding a duplicate. BeginPackage declares all sub-packages as
   dependencies; EndPackage leaves all their contexts on $ContextPath, making
   every public symbol accessible unqualified after Needs["MFGraphs`"].

   Package loading DAG (each ↓ leads from the line directly above it; every
   package also needs primitives):
     primitives          (shared primitives: j, u, z, alpha, Cost, flags)
       ↓
     utilities           (typed-object infra, rule management, shared helpers)
       ↓
     scenarioTools       (typed scenario kernel)
       ↓
     examples            (named scenario factories)
     unknownsTools       (symbolic unknown families)
       ↓
     systemTools         (structural equation system)
       ↓
     validationTools     (exact soundness/completeness checks)
     graphicsTools       (visualization helpers)
     solversTools        (reduceSystem and the *ReduceSystem reducers)
       ↓
     orchestrationTools  (solveScenario: high-level orchestration)

   Optional subpackages (NOT loaded here; load explicitly after Needs["MFGraphs`"]):
     Tawaf`         (unrolled circumambulation scenario builder)
     numericOracle` (LP-relaxation oracle; introduces FindInstance floats)
*)

(* Pure function rather than With[{dir = ...}] so no helper symbol is created
   in the caller's context. *)
$Path = (Prepend[DeleteCases[$Path, #], #] &)[DirectoryName[$InputFileName]];

BeginPackage["MFGraphs`",
  { "primitives`",
    "utilities`",
    "scenarioTools`",
    "examples`",
    "unknownsTools`",
    "systemTools`",
    "validationTools`",
    "solversTools`",
    "orchestrationTools`",
    "graphicsTools`"
  }
];

EndPackage[]
