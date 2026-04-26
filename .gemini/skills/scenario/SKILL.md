---
name: scenario
description: Enforce the canonical MFGraphs scenario-kernel pipeline from scenario construction through ReduceSystem solving.
---

# /scenario

Use this exact workflow:

1. Construct
```mathematica
s = makeScenario[rawAssociation]
```

2. Unknowns
```mathematica
unk = makeUnknowns[s]
```

3. System
```mathematica
sys = makeSystem[s, unk]
```

4. Validation
```mathematica
scenarioQ[s]
mfgSystemQ[sys]
```
Both must be `True` before proceeding.

5. Consistency Check (before solve)
```mathematica
SystemData[sys, "consistentCosts"]
```
Switching-cost consistency must be satisfied before solving.

6. Solve
```mathematica
d2e = SystemDataFlatten[sys]
result = ReduceSystem[sys]
```
