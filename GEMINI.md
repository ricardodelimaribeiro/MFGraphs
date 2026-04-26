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
s = GetExampleScenario[12][{{1, 100}}, {{4, 0}}, {}];

(* Generate symbolic unknowns and structural equations *)
unk = makeUnknowns[s];
sys = makeSystem[s, unk];

(* Access structural data *)
data = SystemDataFlatten[sys];
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
| `MFGraphs/Scenario.wl` | Typed scenario kernel |
| `MFGraphs/Unknowns.wl` | Symbolic unknown bundle construction |
| `MFGraphs/System.wl` | Structural equation system kernel |
| `Scripts/RunTests.wls` | Test suite runner |
