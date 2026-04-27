# Troubleshooting

## Solver and Performance Issues

### Solver timeout during testing

**Symptom:** Test case times out (exceeds `60s`, `300s`, `900s`, or `1800s` depending on tier)

**Causes & Solutions:**
1. **Running a large case** — Check which tier you're in (see `BENCHMARKS.md` § "Test case tiers")
   - Small tier: 60s timeout (cases 1–6, 27)
   - Core tier: 300s timeout (cases 9, 10, 12, 14, 15, 18, Paper)
   - Large tier: 900s timeout (cases 13, 19–23, Braess variants, Grid0303)
   - VLarge tier: 1800s timeout (Grid1020 only — stress test)

2. **Enable verbose logging** to see progress:
   ```mathematica
   $MFGraphsVerbose = True;
   result = CriticalCongestionSolver[d2e];  (* prints timing per stage *)
   ```

3. **Try CriticalCongestionSolver first** (much faster than NonLinearSolver):
   ```mathematica
   critical = CriticalCongestionSolver[d2e];  (* ~seconds *)
   nonlinear = NonLinearSolver[d2e];          (* ~minutes *)
   ```

4. **Check if DNFReduce is the bottleneck** — enables `BottleneckReport.wls`:
   ```bash
   wolframscript -file Scripts/BottleneckReport.wls
   # See Results/bottleneck_report.md for timing breakdown
   ```

### Solver returns "Infeasible"

**Symptom:** `result["Status"]` is `"Infeasible"` or `IsFeasible[result]` returns `False`

**Causes & Solutions:**

1. **Switching costs violate triangle inequality** (expected failure):
   ```mathematica
   IsSwitchingCostConsistent[switchingCosts]  (* should return True *)
   ```
   If False, switching costs form cycles with negative "profit" — the problem has no feasible equilibrium. **This is expected behavior** and indicates the feature is working correctly.

2. **Flow variables contain negative values** — indicates no feasible non-negative flow solution exists:
   ```mathematica
   result = CriticalCongestionSolver[d2e];
   If[!IsFeasible[result],
     Print["Negative flows: ", Select[result["AssoCritical"], # < 0 &]]
   ]
   ```

3. **Stale cache from previous problem** — clear the solver cache:
   ```mathematica
   result1 = CriticalCongestionSolver[d2e1];
   ClearSolveCache[];  (* IMPORTANT: clear between different problems *)
   result2 = CriticalCongestionSolver[d2e2];
   ```
   Failure to clear can cause solutions from one problem to leak into another.

### RecursionLimit exceeded

**Symptom:** Error message contains "RecursionLimit" or computation stops abruptly

**Expected for certain cases:**
- **Grid1020** (vlarge tier) is a deliberate stress test designed to hit `RecursionLimit` by design
- Used to validate solver behavior under recursion constraints

**For other cases** (unexpected RecursionLimit):
1. Check network complexity — very deep Or-branch nesting in DNFReduce can trigger this
2. Increase recursion limit (temporary workaround, not recommended):
   ```mathematica
   $RecursionLimit = 100000;  (* default is 4096 *)
   ```
3. Simplify the network or contact developers for optimization suggestions

## Configuration and Parallel Issues

### Parallel kernels not launching

**Symptom:** Solver runs serially despite list length > 6

**Causes & Solutions:**

1. **Verify parallel threshold is exceeded:**
   ```mathematica
   $MFGraphsParallelThreshold = 6  (* default *)
   (* Your operation needs list length > 6 to trigger parallel dispatch *)
   ```

2. **Ensure kernels are initialized and definitions distributed:**
   ```mathematica
   EnsureParallelKernels[];
   DistributeDefinitions[alpha, g, V];
   $MFGraphsParallelReady = True;
   ```

3. **Enable verbose to confirm dispatch:**
   ```mathematica
   $MFGraphsVerbose = True;
   result = NonLinearSolver[d2e];  (* will print timing, including parallel dispatch *)
   ```

4. **Check for solver parameter overrides** — if you've overridden `alpha`, `g`, or `V`, they must be distributed to all kernels via `DistributeDefinitions`

### Incorrect or stale parallel results

**Symptom:** Parallel results differ from serial results, or results don't reflect latest code changes

**Cause:** `$SolveCache` is not shared across parallel kernels — each kernel has its own cache

**Solution:**
1. Call `ClearSolveCache[]` before switching from serial to parallel solving
2. After making code changes to `DNFReduce`, call `ClearSolveCache[]` and restart parallel kernels:
   ```mathematica
   ClearSolveCache[];
   ParallelEvaluate[ClearSolveCache[]];  (* clear on all kernels *)
   DistributeDefinitions[...];  (* redistribute updated definitions *)
   ```

## Testing Issues

### Test failures with "expected failure" messages

**Symptom:** Test output shows `FAILED` but indicates it's expected

**What this means:** The test validates that certain problems **correctly fail**. Examples:
- Switching costs that violate triangle inequality (problem has no valid equilibrium)
- Network configurations that produce infeasible systems

**Action:** No action needed — these are feature validations, not bugs

### Some tests timeout but others pass

**Symptom:** Fast suite (`RunTests.wls fast`) passes, but slow suite has timeouts

**Solution:** Slow suite includes large-case solvers that take minutes. Check:
1. Is your machine under heavy load?
2. Are you running other Mathematica sessions simultaneously?
3. Can you run a specific slow case with a longer timeout?
   ```bash
   wolframscript -file Scripts/BenchmarkSuite.wls large case=13 timeout=1200
   ```

## Data Format Issues

### "Invalid network data" or unexpected solver errors

**Symptom:** Error mentioning network data keys or association format

**Verify network data format:**
```mathematica
Data = <|
  "Vertices" -> {1, 2, 3, 4},
  "Adjacency" -> {{0,1,0,0}, {0,0,1,0}, {0,0,0,1}, {0,0,0,0}},
  "Entries" -> {{1, 200}},
  "Exits" -> {{4, 0}},
  "Switching" -> {{1,2,3, 5}, {3,2,1, 5}}
|>
```

**Common mistakes:**
- Using `Lists` instead of `Associations` for network data
- Missing a required key
- Mismatched vertex labels (e.g., adjacency matrix 0-indexed but vertices list 1-indexed)
- Entering switching costs in wrong order (should be `{from, at, to, cost}`)

### Symbolic parameters not substituted

**Symptom:** Solver fails or returns symbolic results containing `I1`, `U1`, `S1`, etc.

**Solution:** Examples use symbolic parameters; substitute before solving:
```mathematica
Data = GetExampleData[12];  (* contains I1, U1 *)
Data = Data /. {I1 -> 100, U1 -> 0};  (* substitute before solving *)
d2e = DataToEquations[Data];
result = CriticalCongestionSolver[d2e];
```

Default parameters for benchmarks are in `Scripts/BenchmarkHelpers.wls` (`$DefaultParams`).

## Getting More Help

### Enable comprehensive debugging

```mathematica
$MFGraphsVerbose = True;
result = CriticalCongestionSolver[d2e];
(* Check result structure: *)
Keys[result]  (* shows what data is available *)
result["Status"]  (* see if Feasible or Infeasible *)
```

### Check solver-specific outputs

```mathematica
(* For CriticalCongestionSolver *)
result["AssoCritical"]  (* zero-flow equilibrium *)

(* For NonLinearSolver *)
result["AssoNonCritical"]  (* general congestion solution *)

(* For MonotoneSolver *)
result  (* either solution or Null *)
```

### Review relevant documentation

- **CLAUDE.md** — Developer guide with architecture and debugging workflows
- **BENCHMARKS.md** — Performance profiling methodology and bottleneck details
- **API_REFERENCE.md** — Complete function signatures and options
- **README.md** — Quick start examples and configuration

### Report or discuss issues

Detailed issues can be reported on GitHub with:
1. Minimal reproducible example (data + solver call)
2. Output of `$MFGraphsVerbose = True` run
3. Mathematica version
4. OS and machine specs (optional but helpful)
