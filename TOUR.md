# MFGraphs Tour

A guided reading order for newcomers. Skim the markdown first, then open the workbooks in Mathematica and evaluate cell-by-cell.

## What MFGraphs is

A Wolfram Language paclet for solving stationary **Mean Field Games on networks** in the critical-congestion regime (`α = 1`). A large population of rational agents moves through a graph; mass conservation (Kirchhoff laws) plus complementarity yields a Linear Complementarity Problem that the package builds and reduces symbolically. See `README.md` for the math.

## Read first (~30 min, in this order)

1. **`README.md`** — model overview, Hamiltonian, mass-conservation / complementarity setup.
2. **`CLAUDE.md`** — the most useful map of the repo: load order, the `makeScenario → makeSymbolicUnknowns → makeSystem → solveScenario` pipeline, the typed-wrapper pattern, naming conventions. Reads like an architecture guide.
3. **`API_REFERENCE.md`** — auto-generated public surface. Skim the headings to know what exists; come back as a reference.
4. **`CONTRIBUTING.md`** + **`TROUBLESHOOTING.md`** — skim only, return when you hit a snag.

## Open in Mathematica, in this order

5. **`MFGraphs/Getting started.wl`** — the intended on-ramp notebook. Open in the front end and evaluate cells top-to-bottom. Covers all 13 capabilities of the core scenario-kernel phase: scenario constructors, named factories, `scenarioData` inspection, per-edge Hamiltonian model, symbolic unknowns, `makeSystem`, chain examples, topology and solution visualisation, the augmented (state-space) graph, and the canonical Jamaratv9 example.
6. **`MFGraphs/TawafWorkbook.wl`** — the most polished worked example. Sections 1.1 → 1.9 walk one ring through unrolling, physical-edge coupling (the package's central modelling idea for shared-pavement scenarios), helix view, and the opt-in density extension. Pair with `docs/research/notes/tawaf-model.md` for the modelling rationale.
7. **`MFGraphs/Jamarat.wl`** — multi-entrance / multi-exit example built on the registered `Jamaratv9` scenario. Demonstrates the value-function gradient via `richNetworkPlot[..., UseColorFunction -> True]`.

## Read source in load order (when you're ready to dig)

The flat-context load order from `CLAUDE.md`:

8. `MFGraphs/primitives.wl` → `MFGraphs/utilities.wl` (foundations)
9. `MFGraphs/scenarioTools.wl` (`makeScenario`, validation, topology cache)
10. `MFGraphs/unknownsTools.wl` (the `j`, `u`, `z` symbolic families)
11. `MFGraphs/systemTools.wl` (boundary / flow / complementarity / Hamiltonian builders)
12. `MFGraphs/solversTools.wl` (DNF + active-set reducers)
13. `MFGraphs/orchestrationTools.wl` (`solveScenario` wrapper)
14. `MFGraphs/graphicsTools.wl` (`rawNetworkPlot`, `richNetworkPlot`)
15. `MFGraphs/Tawaf.wl` (the unroll-then-couple specialisation — read after the core)

## Tests as documentation

16. `MFGraphs/Tests/scenario-kernel.mt` and `MFGraphs/Tests/scenario-consistency.mt` are the cleanest worked specs of what each stage guarantees — short, readable, authoritative when the prose drifts.

## The one-liner

Open **`MFGraphs/Getting started.wl`** in Mathematica with **`README.md`** and **`CLAUDE.md`** open in another window. Everything else is reference.
