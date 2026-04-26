---
name: scenario
description: Walk the MFGraphs solver pipeline: makeScenario → makeSystem → ReduceSystem with consistency checks
---

You are executing the `/scenario` workflow for the MFGraphs project. This skill documents the canonical pipeline for constructing and solving a Mean-Field Game on a graph.

## Pipeline

### Step 1 — Build the scenario
```mathematica
s = makeScenario[raw]
```
`raw` is an Association with at minimum:
- `"Model"` → `<| "Vertices List" → ..., "Adjacency Matrix" → ..., "Entrance Vertices and Flows" → ..., "Exit Vertices and Terminal Costs" → ..., "Switching Costs" → ... |>`

Alternatively use a factory: `s = GetExampleScenario[key][entries, exits, sc, alpha, V, g]`

### Step 2 — Validate the scenario
```mathematica
scenarioQ[s]   (* must return True *)
```
If this returns `False`, inspect `ScenarioData[s]` for missing or invalid keys and fix the input before proceeding.

### Step 3 — Build unknowns and system
```mathematica
unk = makeUnknowns[s]
sys = makeSystem[s, unk]
mfgSystemQ[sys]   (* must return True *)
```

### Step 4 — Check cost consistency
```mathematica
d2e = SystemDataFlatten[sys]
d2e["consistentCosts"] === True
```
If `consistentCosts` is not `True`, the switching-cost matrix violates the triangle inequality on some triple. Inspect `d2e` for the offending entries and correct the scenario's `"Switching Costs"` before proceeding.

### Step 5 — Solve
```mathematica
result = ReduceSystem[sys]
```
`ReduceSystem` calls Mathematica's `Reduce` on the full linear/complementarity system. It returns an Association with keys including `"Solutions"`, `"TimingSeconds"`, and `"BoundData"`.

Note: `ReduceSystem` does not incorporate switching costs (no `AltOptCond`). For switching-cost-aware solving, the planned approach is a two-pass strategy: solve without `AltOptCond` first, then reintroduce it.

## Interpreting results
- `result["Solutions"]` — list of solution rules `{j[...] → ..., u[...] → ..., ...}`
- `result["TimingSeconds"]` — wall-clock time
- If the list is empty, the system may be infeasible or `Reduce` timed out; check `$MFGraphsVerbose = True` for diagnostics.
