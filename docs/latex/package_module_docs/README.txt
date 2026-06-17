Per-module TeX documentation for the MFGraphs package.

Active package surface (loaded by Needs["MFGraphs`"]):
  MFGraphs_loader.tex   -- top-level loader
  primitives.tex        -- symbolic atoms and runtime flags
  utilities.tex         -- shared typed-object and rule helpers
  scenarioTools.tex     -- typed scenario kernel
  examples.tex          -- scenario constructors and example registry
  unknownsTools.tex     -- symbolic unknown bundle
  systemTools.tex       -- structural equation system
  solversTools.tex      -- symbolic solvers and validation
  orchestrationTools.tex -- solveScenario, SolveMFG, clearSolveCache
  graphicsTools.tex     -- rawNetworkPlot, richNetworkPlot, augmentAuxiliaryGraph

Opt-in subpackages (NOT loaded by Needs["MFGraphs`"]; load explicitly):
  Tawaf.tex             -- Needs["Tawaf`"]: unrolled circumambulation builder
  numericOracle.tex     -- Needs["numericOracle`"]: LP-relaxation oracle + solveScenarioWithOracle

Build any file with:
  pdflatex <file>.tex

The .tex sources are the long-form companion to API_REFERENCE.md (auto-generated
from ::usage strings) and to CLAUDE.md (architectural overview).
