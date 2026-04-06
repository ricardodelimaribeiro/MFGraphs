# Phase 0: Jamaratv9 Analysis Findings

## 0.1: Edge Count & Adjacency Reconciliation

**Status**: Ready to verify

Jamaratv9 data from ExamplesData.wl (lines 302-318):
```
Vertices: {1,2,3,4,5,6,7,8,9}
Adjacency matrix: 9x9 matrix with entries
```

From manual inspection of adjacency matrix:
```
{0,1,1,0,0,0,0,0,0},  (* 1→2, 1→3 *)
{0,0,0,1,0,0,0,0,0},  (* 2→4 *)
{0,0,0,1,1,0,0,0,0},  (* 3→4, 3→5 *)
{0,0,0,0,0,1,0,0,0},  (* 4→6 *)
{0,0,0,0,0,1,1,0,0},  (* 5→6, 5→7 *)
{0,0,0,0,0,0,0,1,0},  (* 6→8 *)
{0,0,0,0,0,0,0,0,1},  (* 7→9 *)
{0,0,0,0,0,0,0,0,1},  (* 8→9 *)
{0,0,0,0,0,0,0,0,0}   (* no edges from 9 *)
```

**Directed edge count**: 11 edges (1→2, 1→3, 2→4, 3→4, 3→5, 4→6, 5→6, 5→7, 6→8, 7→9, 8→9)

**Paper claim**: "14 edges" (need to check if this includes augmented reverse edges or auxiliary vertices in the expanded graph)

**Next step**: Run GetExampleData["Jamaratv9"] in Mathematica and count edges programmatically; cross-check with `Jama.pdf` figure

---

## 0.2: Extract True `q` from Preprocessing

**Status**: Requires running DataToEquations in Mathematica

Plan:
```mathematica
Get["MFGraphs/MFGraphs.wl"];
Get["Scripts/BenchmarkHelpers.wls"];
data = GetExampleData["Jamaratv9"] /. GetParams["Jamaratv9"];
d2e = DataToEquations[data];
(* Print InitRules to see fixed sign assignments *)
Print["InitRules: ", d2e["InitRules"]];
(* Print NewSystem to count Or-branches (uncertain triples) *)
Print["NewSystem: ", d2e["NewSystem"]];
(* Expected: verify 2^95 or measure actual q *)
```

---

## 0.3: Function Attributes

**Status**: Verified (no special attributes)

From code inspection of DataToEquations.wl and DNFReduce.wl:
- `CriticalCongestionSolver[Eqs_]` — no HoldFirst/HoldAll
- `TripleClean[{{EE_, NN_, OR_}, rules_}]` — no HoldFirst/HoldAll
- `TripleStep[{{EEs_, NNs_, ORs_}, rules_List}]` — no HoldFirst/HoldAll
- `DNFReduce[xp_, sys_]` — multiple overloads, no HoldFirst/HoldAll

**Wrapping strategy**: Safe to use `args___` splat pattern

---

## 0.4: Pruning Sites in DNFReduce.wl

**Status**: Annotated

### Immediate pruning (False branches)
- **Location**: DNFReduce.wl lines 101-109 (While loop in conjunction handling)
- **Trigger**: `If[xpAcc === False, Return[False, Module]]` and `If[elem === False, Return[False, Module]]`
- **Counter**: `branchesPrunedImmediate++` when False is returned immediately

### Or-branch short-circuit
- **Location**: DNFReduce.wl lines 157-172 (Or handling)
- **Trigger**: `If[result1 === xp, ...]` — if first Or-branch doesn't reduce the space, we exit early
- **Counter**: This is an **optimization**, not a pruning event per se (branch already counted as "entered")

### Subsumption pruning (implication checks)
- **Location**: DNFReduce.wl lines 222-256 (SubsumptionPrune)
- **Trigger**: Pairwise `Implies[branches[[i]], branches[[j]]]` checks
- **Counter**: `branchesPrunedElimination++` when `keep[[j]] = False`
- **Called from**: ReduceDisjuncts (line 274, 280)

### Post-reduction via ReduceDisjuncts
- **Location**: DNFReduce.wl lines 262-299
- **Trigger**: After DNFReduce returns a DNF, ReduceDisjuncts applies Reduce or SubsumptionPrune
- **Counter**: Track input vs output branch count; difference is elimination

---

## 0.5: Parameter Verification

**Status**: Verified from BenchmarkHelpers.wls

From `Scripts/BenchmarkHelpers.wls` lines 123-162:

Jamaratv9 is **not** in `$CaseParams`, so it uses `$DefaultParams`:

```mathematica
$DefaultParams = {
  I1 -> 100,    (* entrance flow at vertex 1 *)
  I2 -> 50,     (* entrance flow at vertex 2 *)
  I3 -> 30,
  U1 -> 0,      (* terminal cost at exit vertex 7 *)
  U2 -> 0,      (* terminal cost at exit vertex 8 *)
  U3 -> 0,      (* terminal cost at exit vertex 9 *)
  S1 -> 1, ..., S16 -> 1  (* switching costs, not used for Jamaratv9 *)
}
```

**Canonical substitution for Jamaratv9**: `{I1 -> 100, I2 -> 50, U1 -> 0, U2 -> 0, U3 -> 0}`

---

## 0.6: Define "Surviving Leaves M"

**Status**: Decision pending

Three possible definitions:

1. **M1**: Branches not eliminated by immediate constraints
   - Measured as: `branchesEntered - branchesPrunedImmediate`
   - Semantic: "Made it through linear feasibility check"

2. **M2**: Algebraic systems that didn't become inconsistent during TripleClean/DNFReduce
   - Measured as: branches that reach ReduceDisjuncts output without being False
   - Semantic: "Algebraically consistent branches"

3. **M3**: Branches yielding unique numerical solutions from NonLinearSolver
   - Measured as: length(result["AssoNonCritical"]) or similar
   - Semantic: "Actually-solvable branches"

**Recommendation**: Instrument to track all three, then footnote in paper which one we report in abstract/conclusion.

---

## Summary Table

| Task | Status | Finding |
|------|--------|---------|
| 0.1: Edge count | Pending | 11 directed edges (verify against "14" claim) |
| 0.2: q value | Pending | Need Mathematica run |
| 0.3: Function attributes | **✓ Verified** | No special attributes; args___ wrapping safe |
| 0.4: Pruning sites | **✓ Annotated** | DNFReduce.wl lines 101-109, 222-256; ReduceDisjuncts lines 262-299 |
| 0.5: Parameters | **✓ Verified** | I1→100, I2→50, U1→0, U2→0, U3→0 |
| 0.6: Surviving leaves def | Pending | Choose between M1, M2, M3; recommend track all three |

---

## Next Steps

1. **Run in Mathematica** (on your machine, not sandbox):
   ```mathematica
   Get["MFGraphs/MFGraphs.wl"];
   Get["Scripts/BenchmarkHelpers.wls"];
   data = GetExampleData["Jamaratv9"] /. GetParams["Jamaratv9"];
   d2e = DataToEquations[data];
   Print["Edge count: ", Total[Flatten[data["Adjacency Matrix"]]]];
   Print["InitRules keys: ", Keys[d2e["InitRules"]]];
   Print["NewSystem: ", d2e["NewSystem"]];
   ```

2. **Cross-check paper claim**: Open `figures/Jama.pdf` and verify edge labels against extracted adjacency

3. **Make M definition choice** and document in instrumentation script

4. Proceed to Phase 1 (parameter table extraction)

