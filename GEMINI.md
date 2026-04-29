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

(* Generate symbolic unknowns and structural equations *)
unk = makeUnknowns[s];
sys = makeSystem[s, unk];

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
| `MFGraphs/scenarioTools.wl` | Typed scenario kernel |
| `MFGraphs/examples.wl` | Example scenario factories |
| `MFGraphs/unknownsTools.wl` | Symbolic unknown bundle construction |
| `MFGraphs/systemTools.wl` | Structural equation system kernel |
| `MFGraphs/solversTools.wl` | Critical-congestion structural solver |
| `MFGraphs/graphicsTools.wl` | Scenario and solution plotting helpers |
| `Scripts/RunTests.wls` | Test suite runner |
