# GEMINI.md - MFGraphs Project Guidance

For comprehensive guidance on working with MFGraphs, see **[CLAUDE.md](CLAUDE.md)**.

## Quick Reference

**MFGraphs** is currently in a **core scenario-kernel phase**.
The active package surface focuses on typed scenario, unknown, and structural system construction.

### Load the package
```mathematica
Needs["MFGraphs`"]
```

### Quick start (Core Scenario Kernels)
```mathematica
(* Build a scenario using a factory or makeScenario *)
s = getExampleScenario[12, {{1, 100}}, {{4, 0}}];

(* Generate exact symbolic unknowns and structural equations *)
unk = makeSymbolicUnknowns[s];
sys = makeSystem[s, unk];

(* Solve through the default DNF-first orchestration path *)
sol = solveScenario[s];

(* Opt in to exact branch-state solvers when comparing solver performance *)
optSol = solveScenario[s, optimizedDNFReduceSystem];
activeSol = solveScenario[s, activeSetReduceSystem];

(* Access structural data *)
data = systemDataFlatten[sys];
data["AltFlows"]
```

## For Detailed Documentation

→ See **[CLAUDE.md](CLAUDE.md)** for:
- Environment setup and prerequisites
- Complete Quick Start with code examples
- Architecture overview (Scenario, Unknowns, System)

→ See **[CONTRIBUTING.md](CONTRIBUTING.md)** for:
- PR workflow and code style
- Testing procedures
- Commit message conventions

## Key Files

| File | Purpose |
|------|---------|
| `MFGraphs/MFGraphs.wl` | Package loader |
| `MFGraphs/primitives.wl` | Symbolic atoms (`j`, `u`, `z`, `alpha`, `Cost`) and runtime flags |
| `MFGraphs/utilities.wl` | Shared typed-object, rule, and critical-congestion helpers |
| `MFGraphs/scenarioTools.wl` | Typed scenario kernel |
| `MFGraphs/examples.wl` | Example scenario factories |
| `MFGraphs/unknownsTools.wl` | Symbolic unknown bundle construction |
| `MFGraphs/systemTools.wl` | Structural equation system kernel |
| `MFGraphs/solversTools.wl` | Critical-congestion structural solvers |
| `MFGraphs/orchestrationTools.wl` | High-level DNF-first solver orchestration |
| `MFGraphs/graphicsTools.wl` | Scenario and solution plotting helpers |
| `MFGraphs/Tawaf.wl` | **Opt-in** unrolled circumambulation scenario builder |
| `MFGraphs/numericOracle.wl` | **Opt-in** LP-relaxation oracle + `solveScenarioWithOracle` |
| `Scripts/RunTests.wls` | Test suite runner |

Opt-in subpackages are not loaded by `Needs["MFGraphs`"]`. Load them
explicitly: `Needs["Tawaf`"]`, `Needs["numericOracle`"]`.

Benchmark selector names include `dnf` (default), `optimizeddnf`, `activeset`,
`reduce`, `boolean`, and `findinstance`.

Result-kind classification is flow-first: primary classification depends on
edge flows `j[_,_]` and values `u[__]`; transition-flow determinacy for
`j[_,_,_]` is reported separately. Both helpers (`solutionResultKind`,
`solutionVariableDiagnostics`) are internal — they live in
`` solversTools`Private` `` and are not part of the exported API; use the
public `solutionReport[sys, sol]` for these diagnostics.

Use `solutionReport[sys, sol]` for read-only solution diagnostics and
`solutionBranchCostReport[sys, sol]` to rank residual branches. The direct and
flow-first critical solvers are explicit opt-ins (`directCriticalSystem`,
`flowFirstCriticalSystem`); the default `solveScenario` path remains
`dnfReduceSystem`.

Hamiltonian signed-cost convention: for physical edge `{a,b}`, use oriented net
flow `m = j[a,b] - j[b,a]`; `m Sign[m]` is the critical traversal cost `Abs[m]`.
Active critical solvers consume `EqGeneral`; `Nrhs`/`CostCurrents` are metadata.
