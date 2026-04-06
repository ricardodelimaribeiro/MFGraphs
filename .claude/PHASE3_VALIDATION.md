# Phase 3: Validation Run (fork/merge test cases)

## Goal
Verify that the instrumentation in `Scripts/InstrumentJamarat.wls` works correctly and doesn't significantly impact performance, before running the 12-minute Jamaratv9 solver.

## Prerequisite
Have Mathematica 12.0+ installed and `wolframscript` in your PATH.

## Test Cases
We'll validate on two small, well-understood cases:
1. **fork1-3**: Expected to have 3 branches, tests branching behavior
2. **merge1-3**: Expected to have 3 branches, tests merging behavior

Both should complete in seconds, so any issues will show up quickly.

## Step 1: Create a minimal validation script

Save this as `Scripts/ValidateInstrumentation.wls`:

```mathematica
#!/usr/bin/env wolframscript
(* Minimal validation of instrumentation on fork/merge examples *)

Get["MFGraphs/MFGraphs.wl"];
Get["Scripts/BenchmarkHelpers.wls"];

(* Load instrumentation wrappers *)
Get["Scripts/InstrumentJamarat.wls"];  (* This redefines CriticalCongestionSolver, etc. *)

(* Test 1: fork1 *)
Print["\n=== VALIDATION TEST 1: fork1 ==="];
prunStats = <|
  "timestamp" -> Now, "network" -> "fork1",
  "q" -> Null, "branchesEntered" -> 0, "branchesPrunedImmediate" -> 0,
  "branchesPrunedElimination" -> 0, "leavesSurviving" -> 0,
  "totalWallTime" -> 0.0, "passTimings" -> {},
  "dnfReduceStats" -> <|"branchesIn" -> 0, "branchesOut" -> 0|>,
  "reduceDisjunctsStats" -> <|"input" -> 0, "output" -> 0|>
|>;

data1 = GetExampleData[7];  (* fork example *)
params1 = GetParams[7];
data1 = data1 /. params1;
d2e1 = DataToEquations[data1];
result1 = CriticalCongestionSolver[d2e1];

Print["\nfork1 Results:"];
Print["  Feasible: ", IsFeasible[result1]];
Print["  Branches entered: ", prunStats["branchesEntered"]];
Print["  Branches survived: ", prunStats["leavesSurviving"]];
Print["  Total time: ", N[prunStats["totalWallTime"], 5], " s"];
Print["  TripleClean passes: ", Length[prunStats["passTimings"]]];

(* Test 2: merge1 *)
Print["\n=== VALIDATION TEST 2: merge1 ==="];
prunStats = <|
  "timestamp" -> Now, "network" -> "merge1",
  "q" -> Null, "branchesEntered" -> 0, "branchesPrunedImmediate" -> 0,
  "branchesPrunedElimination" -> 0, "leavesSurviving" -> 0,
  "totalWallTime" -> 0.0, "passTimings" -> {},
  "dnfReduceStats" -> <|"branchesIn" -> 0, "branchesOut" -> 0|>,
  "reduceDisjunctsStats" -> <|"input" -> 0, "output" -> 0|>
|>;

data2 = GetExampleData[9];  (* merge example *)
params2 = GetParams[9];
data2 = data2 /. params2;
d2e2 = DataToEquations[data2];
result2 = CriticalCongestionSolver[d2e2];

Print["\nmerge1 Results:"];
Print["  Feasible: ", IsFeasible[result2]];
Print["  Branches entered: ", prunStats["branchesEntered"]];
Print["  Branches survived: ", prunStats["leavesSurviving"]];
Print["  Total time: ", N[prunStats["totalWallTime"], 5], " s"];
Print["  TripleClean passes: ", Length[prunStats["passTimings"]]];

Print["\n=== VALIDATION COMPLETE ==="];
```

## Step 2: Run the validation script

From the repo root in your Mathematica environment:

```bash
wolframscript -file Scripts/ValidateInstrumentation.wls
```

Expected output: Both tests should complete in < 5 seconds each. Check that:
- ✓ No errors or warnings
- ✓ Both show `Feasible: True` (or expected feasibility)
- ✓ `branchesEntered` and `leavesSurviving` are reasonable (small for these simple cases)
- ✓ `TripleClean passes` are non-zero (≥ 1)
- ✓ Total time < 10 seconds for both

## Step 3: Verify instrumentation overhead

Compare with an uninstrumented run:

```mathematica
(* In Mathematica, without instrumentation *)
Get["MFGraphs/MFGraphs.wl"];
Get["Scripts/BenchmarkHelpers.wls"];

data = GetExampleData[7] /. GetParams[7];
d2e = DataToEquations[data];
{t1, result1} = AbsoluteTiming[CriticalCongestionSolver[d2e]];
Print["Without instrumentation: ", t1, " s"];
```

Then compare with the instrumented time from Step 2. Overhead should be < 10%.

## Troubleshooting

If the validation script fails:

1. **Error: symbol not defined** → Make sure `Get["Scripts/InstrumentJamarat.wls"]` is executed before solver calls
2. **Error: Division by zero or unexpected Null** → The wrapper may be interfering with return values; check the wrapping logic
3. **Significant slowdown (> 20%)** → Switch from `AppendTo` to `Sow`/`Reap` in InstrumentJamarat.wls
4. **Counters are zero** → The wrapper functions may not be evaluating; check that they're actually replacing the originals

## Proceed to Phase 4 only after:
- [ ] ValidateInstrumentation.wls runs without errors
- [ ] fork1 and merge1 show expected branching behavior
- [ ] Overhead < 10%
- [ ] All counters (branchesEntered, leavesSurviving, passTimings) are non-zero and reasonable
