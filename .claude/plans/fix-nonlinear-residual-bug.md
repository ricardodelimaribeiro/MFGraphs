# Plan: Fix the silent non-linear residual bug in NonLinearSolver and MonotoneSolver

## TL;DR

Both `NonLinearSolver` and `MonotoneSolver` currently report success on every α-non-trivial test case, but their solutions do **not** satisfy the actual non-linear MFG equations `Nlhs == Nrhs`. The reported `KirchhoffResidual` only checks linear flow balance (`KM·j = B`); the equation residual `|Nlhs - Nrhs|∞` is 25–50 across all Y7 / Case12 / Braess test cases — i.e. **the equation is grossly violated**, but no surface in the package detects it.

This plan fixes the bug in **six phases**, ordered to (a) build detection before fix, (b) isolate root cause before patching, and (c) lock in regression coverage before declaring victory.

---

## What we know (evidence base)

1. **Symptom** (`Scripts/VerifyNonLinearSolutions.wls` runs):
   | Case | α | Kirchhoff | `|Nlhs−Nrhs|∞` |
   |------|---|-----------|----------------|
   | Y7_alpha2 | 2   | 0 | **50.0** |
   | Y7_alpha125 | 1.25 | 0 | **50.0** |
   | Y7_alpha15 | 1.5 | 0 | **50.0** |
   | Case12_alpha25 | 2.5 | 0 | **25.0** |
   | Case12_alpha3 | 3 | 0 | **25.0** |

2. **Direct observation** (`Scripts/DebugNonLinearResidual.wls`): for Y7 with `I1=50`, the solver returns `j[1,2]=50, j[2,1]=0` (Kirchhoff-balanced) but **all 12 `u[…]` variables = 0**. The equation reduces to `u[1,2] - u[2,1] = -50`, so `u` should *not* be zero.

3. **Code-level root cause hypothesis** — there are **two** independent bugs that compose:

   - **Bug A — preprocessing pins `u` to 0** (`MFGraphs/DataToEquations.wl:455-470`): when all switching costs are zero, `MFGPreprocessing` collapses `IneqSwitchingByVertex` to equalities `u[a,b] == u[c,b]`. Combined with `RuleExitValues` (which sets `u[*, exit] = U_exit`) and `EqValueAuxiliaryEdges`, `TripleClean` propagates `u = 0` through the **entire graph** before `MFGSystemSolver` ever runs. `EqGeneral` is then powerless to constrain `u` because `InitRules` already fixes it.

   - **Bug B — convergence test ignores the equation residual** (`MFGraphs/NonLinearSolver.wl:130-136`): the `NestWhileList` predicate is `Norm[Δj]∞ > tol`, i.e. only the change in `j` between iterations. Since the seed comes from `CriticalCongestionSolver` with already-correct `j` values, `flowDelta = 0` after iteration 1, the loop exits as `"ToleranceMet"`, and the wrong `u` is never corrected. The internally printed `"Max error for non-linear solution"` (line 222) is logged via `MFGPrint` but never inspected.

   - **Bug C — `MonotoneSolver` solution is missing `u` and aggregated `j` variables**: `Scripts/VerifyNonLinearSolutions.wls` returns `Missing[NoNumeric]` for every Monotone case because `AssoMonotone` only contains the Kirchhoff-edge flows used by the HRF reduction, not the auxiliary `u[*]`/`j[*,*,*]` variables that `Nlhs`/`Nrhs` reference. Even if Monotone's `j` were correct, the residual cannot be evaluated against its solution dictionary.

   - **Bug D — `BuildSolverComparisonData` reports only the linear residual** (`MFGraphs/MFGraphs.wl:117-122`): the comparison data exposed to every solver result has `KirchhoffResidual` but no `NonLinearResidual`. Every benchmark, every CSV report, and every `IsFeasible` check has been blind to the actual MFG equation.

   - **Bug E — `IsNonLinearSolution` does not gate on the residual** (`MFGraphs/NonLinearSolver.wl:228-269`): it `Reap`s the assoc, prints a "Max error" diagnostic to the screen, and returns the assoc unconditionally. The function name implies a boolean check, but it's actually a verbose printer.

4. **Test coverage gap**: `MFGraphs/Tests/solver-contracts.mt` tests `NonLinearSolver` and `MonotoneSolver` only with `"CongestionExponentFunction" -> Function[edge, 1]` — i.e. the linear case where the bug **happens to be invisible** because at α=1 the equation `u[a,b] - u[b,a] = -j` happens to also be satisfied by certain trivial solutions. There is currently **no automated test that exercises α ≠ 1**, which is why this bug survived.

## What we don't yet know (assumptions to verify in Phase 1)

- Whether the bug also affects α=1 (the linear baseline). The fact that `IntegratedMass[50, edge] ≈ 50` in the verification output suggests `Nrhs ≈ 0` regardless of α — meaning the *same* `u`-zeroing pathology might already corrupt the linear case, and only the absence of an equation-residual check is hiding it.
- Whether `MFGPreprocessing`'s `u`-equality propagation is *correct* for α≠1 (it might be, since `u[a,b] = u[c,b]` is a true value-function statement) and only the *value 0* propagation is wrong.
- Whether `MFGSystemSolver` is *capable* of computing `u` from `j` if `InitRules` did not already pin `u`. (My reading says yes — `EqGeneral /. Ncpc` plus `TripleClean` should solve for `u` given `j` — but this needs empirical confirmation.)
- Whether `MonotoneSolver` has internal access to the value functions (`u`) at the end of the HRF flow, or whether they have to be reconstructed from `j` by post-hoc inversion of the same equation.

These uncertainties shape Phase 1: **gather evidence before fixing anything**.

---

## Phase 0 — Snapshot and branch (10 min)

**Goal**: capture the current broken state in a way we can diff against, then create a branch.

**Actions**:
1. Capture the current verification output to `Results/baseline_broken_residuals.md` with a clear note explaining what each row means.
2. Create a working branch `fix/nonlinear-residual-bug` (current worktree is `vibrant-aryabhata`; branch off whatever its base is).
3. Commit the existing scratch artifacts (`Scripts/VerifyNonLinearSolutions.wls`, `Scripts/DebugNonLinearResidual.wls`, the broken benchmark Markdown files) so the bug is *recorded in git history*, not just in memory.

**Exit criteria**: branch exists, baseline residual table is committed, `git log` shows a commit titled `chore: capture non-linear residual bug evidence`.

**Risk**: none — purely additive.

---

## Phase 1 — Detection infrastructure (must come before any fix) (~2 h)

**Goal**: every solver result and every test must be able to *see* the non-linear residual. After Phase 1 the bug becomes loud and obvious; only then do we earn the right to fix it.

### 1.1  Add `NonLinearResidual` to `BuildSolverComparisonData`

**File**: `MFGraphs/MFGraphs.wl` (lines 86-132)

**Change**: extend the returned association with two new keys:

```mathematica
"NonLinearResidual" -> nonLinResidual,        (* Max[Abs[Nlhs - Nrhs]] *)
"NonLinearResidualVector" -> nonLinDiffs      (* per-equation residuals for diagnostics *)
```

The computation reuses the helper from `Scripts/VerifyNonLinearSolutions.wls`:

```mathematica
ComputeNonLinearResidual[Eqs, solution] := Module[{Nlhs, Nrhs, lhs, rhs, diffs},
    Nlhs = Lookup[Eqs, "Nlhs", $Failed];
    Nrhs = Lookup[Eqs, "Nrhs", $Failed];
    If[Nlhs === $Failed || Nrhs === $Failed, Return[{missing, missing}]];
    lhs = Quiet @ Check[N[Nlhs /. Normal[solution]], $Failed];
    rhs = Quiet @ Check[N[Nrhs /. Normal[solution]], $Failed];
    If[lhs === $Failed || rhs === $Failed, Return[{missing, missing}]];
    diffs = Flatten[{lhs - rhs}];
    diffs = Select[diffs, NumericQ];
    If[diffs === {}, Return[{missing, missing}]];
    {Max[Abs[diffs]], diffs}
];
```

**Important**: the function must accept *any* solution dict shape (the `Normal[]` wrap handles both `Association` and `{_Rule..}`). It must also handle the `MonotoneSolver` case where `solution` doesn't contain `u[*]` keys — return `Missing["MissingVariables"]` (not a numeric number) so consumers can distinguish "violated" from "uncheckable".

### 1.2  Surface the new field in every solver's result

**Files**:
- `MFGraphs/DataToEquations.wl` — `CriticalCongestionSolver` calls `BuildSolverComparisonData` already (line 659, 704); the new field will appear automatically.
- `MFGraphs/NonLinearSolver.wl` — same (line 182).
- `MFGraphs/Monotone.wl` — same (multiple call sites starting line 478).

No code change needed beyond Phase 1.1, as long as `BuildSolverComparisonData` is the single source of truth.

### 1.3  Make `IsFeasible` aware of the non-linear residual

**File**: `MFGraphs/MFGraphs.wl` (find `IsFeasible[result_]` definition)

**Change**: add an *optional* tolerance parameter and treat a numeric `NonLinearResidual` larger than the tolerance as `Infeasible`:

```mathematica
$DefaultNonLinearTolerance = 10^-4;

IsFeasible[result_Association, opts:OptionsPattern[]] :=
  Module[{statusOk, residual, tol},
    statusOk = (* existing logic on Status / Feasibility *);
    residual = Lookup[result, "NonLinearResidual", Missing["NotAvailable"]];
    tol = OptionValue[{Tolerance -> $DefaultNonLinearTolerance}, Tolerance];
    statusOk && (
      !NumericQ[residual] ||      (* uncheckable: trust the legacy status *)
      residual <= tol
    )
  ];
```

**Compatibility**: existing call sites that don't pass `Tolerance` still work; legacy `"Status" -> "Feasible"` results without a numeric `NonLinearResidual` (e.g. Monotone's missing-variables case) keep their old behavior. **Only** results that *can* be checked and *fail* the check flip from `True` to `False`.

### 1.4  Replace `IsNonLinearSolution` with a real predicate

**File**: `MFGraphs/NonLinearSolver.wl` (lines 228-269)

**Change**: split into two:
- `ShowNonLinearDiagnostics[result_]` — keep the existing print-everything behavior under a clearly diagnostic name.
- `IsNonLinearSolution[result_, tol_:10^-4]` — new function that returns `True` only when **all** structural restrictions are `True` *and* `NonLinearResidual <= tol`.

Update the existing usage docstring on line 228 (currently the same name returns the assoc — keep backward compat by having the new predicate return `True/False` and a new diagnostic function carry the old behavior).

### 1.5  Add a regression test that *fails on master*

**File**: `MFGraphs/Tests/nonlinear-residual.mt` (NEW)

```mathematica
(* Regression: NonLinearSolver and MonotoneSolver must satisfy the
   actual non-linear MFG equation, not just flow balance. *)

Test[
    Module[{data, d2e, result},
        data = GetExampleData[7] /. {I1 -> 50, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        result = Quiet @ NonLinearSolver[d2e,
            "CongestionExponentFunction" -> Function[edge, 2],
            "MaxIterations" -> 30, "Tolerance" -> 10^-8];
        NumericQ[result["NonLinearResidual"]] &&
        result["NonLinearResidual"] < 10^-4
    ],
    True,
    TestID -> "Y7 alpha=2 satisfies non-linear MFG equation"
]

(* Same test for α=1.5, α=2.5, plus Case12 alpha=2.5. *)
(* Repeat for MonotoneSolver — *expected to FAIL initially* until Phase 4. *)
```

**Crucial**: this test *must* fail on master (Phase 0 commit) and *must* pass after Phase 3 (NonLinearSolver) and Phase 4 (MonotoneSolver). The diff in test results between phases is the proof we actually fixed something.

### 1.6  Add the new field to the benchmark report

**File**: `Scripts/NonLinearBenchmark.wls`

**Change**: extend `caseResults` and the markdown table to include `NonLinearResidual` as a top-level column. The Phase 0 commit captured the broken Markdown; this commit gives us a *truthful* benchmark that exposes the bug going forward.

**Phase 1 exit criteria**:
- `solver-contracts.mt` still passes (additive change, no API break).
- `nonlinear-residual.mt` exists and **fails** on every NonLinearSolver case with the clear error `NonLinearResidual ≈ 50, expected < 1e-4`.
- A fresh `wolframscript -file Scripts/NonLinearBenchmark.wls` produces a Markdown report whose `NonLinearResidual` column is *non-zero* on every case.
- `IsFeasible` returns `False` on at least one previously-"feasible" verification result.

**Risk**: an `IsFeasible` change could cascade through downstream callers. Mitigation: grep for `IsFeasible[` across the repo before committing; for any caller that depends on the legacy flow-balance-only meaning, add a `Tolerance -> Infinity` override to preserve current behavior. Document explicitly which call sites are intentionally being made stricter.

---

## Phase 2 — Isolate the root cause empirically (~1 h)

**Goal**: between the two bug hypotheses (A: preprocessing pins `u`, B: convergence ignores residual) determine *which* fix is necessary, *which* is sufficient, and *whether* there is a third unknown bug.

**Diagnostic experiments** (each is a small `Scripts/DiagnoseStep*.wls`):

1. **Diag-1**: After `MFGPreprocessing[d2e]`, inspect `InitRules` directly. List every `u[*]` key and its value. Hypothesis A predicts every `u` is `0`.

2. **Diag-2**: Build a *fresh* `d2e` whose preprocessing skips the zero-switching-cost shortcut (force the `Else` branch on `MFGraphs/DataToEquations.wl:471-480`). Run `MFGSystemSolver` with that `d2e` and inspect the resulting `u` values. Two possible outcomes:
   - `u` is now non-zero and satisfies the equation → Bug A is the *only* root cause.
   - `u` is still zero → there is a deeper bug in `MFGSystemSolver`'s handling of `EqGeneral` substitution; need to dig further.

3. **Diag-3**: Take the broken solution (`u=0`, `j=correct`), substitute it into `EqGeneral /. Ncpc` (where `Ncpc` is computed from the actual `j`), and check whether the resulting equation is satisfied. This isolates whether the *system* is being violated or only our *interpretation* of it is.

4. **Diag-4**: Run the α=1 baseline through the same `NonLinearResidual` measurement (now visible thanks to Phase 1). If the linear baseline *also* shows residual ≈ 50, then this is a much wider bug that has been mis-reporting feasibility for the entire history of the project — and a much more careful fix is needed.

5. **Diag-5**: For `MonotoneSolver`, post-process its `j` solution by reconstructing `u` via the equation `u[a,b] - u[b,a] = -sign(Δj)·IntegratedMass(Δj, edge)` and check whether the *enriched* solution satisfies the equation. If yes, Monotone's flow is correct — the only problem is that it doesn't expose `u`.

**Phase 2 exit criteria**: a single Markdown note `Results/diagnosis_nonlinear.md` that names the root cause (A only / B only / A+B / something else) and recommends the precise file edits for Phase 3 and 4. **No code fixes yet** — this phase is reconnaissance.

**Risk**: it is possible that Diag-2 reveals `MFGSystemSolver` cannot solve for `u` even with `InitRules` empty, because `EqGeneral` only constrains *differences* `u[a,b] - u[b,a]` and the system is rank-deficient by 1 per connected component. If so, the fix is more substantial — Phase 3 would need to add a *gauge condition* (e.g. `u[exit, *] = U_exit`) to anchor the solution. The diagnostic should explicitly look for this.

---

## Phase 3 — Fix `NonLinearSolver` (~3 h)

**Precondition**: Phase 2 confirmed root cause(s).

### 3.1  Convergence test (Bug B fix — definitely needed)

**File**: `MFGraphs/NonLinearSolver.wl:127-141`

**Change**: replace the `flowDelta`-only convergence predicate with a residual-aware one:

```mathematica
NestWhileList[
    nonLinearStep[PreEqs],
    AssoNonCritical,
    Function[{prev, curr},
        Module[{flowDelta, eqResidual},
            flowDelta = Norm[N[Values[KeyTake[curr, js]] - Values[KeyTake[prev, js]]], Infinity];
            eqResidual = ComputeNonLinearResidual[PreEqs, curr];
            (* Continue iterating if EITHER measure is above tolerance *)
            (flowDelta > tol) ||
            (NumericQ[eqResidual] && eqResidual > tol)
        ]
    ],
    2, MaxIter
]
```

The `stopReason` logic must be updated to distinguish:
- `"FlowDeltaMet"` — old behavior, but only valid when residual is also small.
- `"ResidualMet"` — both measures small, real success.
- `"NonConverged"` — flow stable but residual still high (the bug we're surfacing).
- `"MaxIterationsReached"` — neither converged.

`resultKind` should be `"Success"` only on `"ResidualMet"`.

### 3.2  Preprocessing fix (Bug A fix — likely needed)

**File**: `MFGraphs/DataToEquations.wl:455-470`

**Hypothesis-dependent change**:

- **If Phase 2 Diag-2 shows that disabling the shortcut fixes the bug**: gate the shortcut on a new condition. The shortcut is sound *only* when the network is also being solved at the linear (α=1) critical case. For non-linear runs we need to keep `IneqSwitchingByVertex` as inequalities and let `MFGSystemSolver` resolve them per-iteration.

  Concretely, the shortcut should not be applied during preprocessing for `NonLinearSolver`. Either:
  - **Option A**: Add a flag `Eqs["__skipUEqualityShortcut"]` that `NonLinearSolver` sets before calling `MFGPreprocessing`.
  - **Option B (preferred)**: Move the equality-extraction logic out of `MFGPreprocessing` and into `DirectCriticalSolver` only (where it is genuinely sound). `MFGPreprocessing` would always keep the inequalities, and `MFGSystemSolver` (which substitutes `Ncpc` per iteration) would resolve them with the actual cost values.

  Option B is the more invasive change but it eliminates the bug at its source rather than papering over it with a flag.

- **If Diag-2 shows the shortcut is fine but TripleClean over-propagates**: the fix is in `TripleClean` itself or in the order of equation application. This is more delicate; mitigation is to leave `u` variables out of `InitRules` until after `Ncpc` substitution.

### 3.3  Gauge anchoring (only if Phase 2 reveals rank deficiency)

If Diag-2 showed that `MFGSystemSolver` produces an under-determined system, add an explicit gauge: pick any vertex `v0` (e.g. an exit) and add `u[v0, *] = U_exit_terminal_cost` to `InitRules` *before* invoking `DNFSolveStep`. This anchors the value function; without it the equation has a one-parameter family of solutions per connected component.

### 3.4  Iteration count and tolerance tuning

After the convergence and preprocessing fixes, the iteration may need more steps to converge for high-α cases. Empirically validate by running:

```bash
wolframscript -file Scripts/NonLinearBenchmark.wls
```

and confirm that all 7 test cases reach `NonLinearResidual < 10^-4` within `MaxIterations -> 50`. If any case fails to converge, investigate the per-iteration residual trajectory (already exposed by `Convergence -> "FinalResidual"`).

**Phase 3 exit criteria**:
- `nonlinear-residual.mt` test cases for `NonLinearSolver` pass.
- `Scripts/VerifyNonLinearSolutions.wls` shows `NL_NonLinear < 10^-4` for every Y7 / Case12 case.
- `Scripts/NonLinearBenchmark.wls` markdown shows non-trivial `NonLinearResidual` columns *all below tolerance*.
- `solver-contracts.mt` still passes (linear case still works — no regression).
- `Scripts/RunTests.wls fast` passes.

---

## Phase 4 — Fix `MonotoneSolver` (~2 h)

**Goal**: `MonotoneSolver` should either (a) emit `u` variables in its solution dict so the residual is computable, or (b) explicitly document and surface the limitation.

### 4.1  Decide: post-hoc reconstruction vs structural fix

The Hessian Riemannian Flow optimizes over `j` directly and never represents `u`. Two approaches:

- **4a (recommended)**: After the HRF flow converges to optimal `j`, reconstruct `u` post-hoc by integrating the equation `u[a,b] - u[b,a] = -sign(Δj)·IntegratedMass(Δj, edge)` along a spanning tree from the exit (where `u = U_exit_terminal_cost`) inward. This gives the *correct* value function corresponding to the equilibrium `j` and is mathematically equivalent to Lagrange multipliers of the flow constraint.

- **4b**: Augment the HRF state to track `u` directly through the flow. More invasive, more correct for time-dependent extensions, but a much larger refactor.

Choose 4a for this PR; 4b is a future enhancement.

### 4.2  Implement `ReconstructValueFunction[d2e_, jSolution_]`

**File**: `MFGraphs/Monotone.wl` (new helper, near `MonotoneVariableFieldValue` line 513)

```mathematica
ReconstructValueFunction[Eqs_Association, jSolution_Association] :=
  Module[{exitVerts, exitCosts, graph, uRules, queue, visited, ...},
    (* Initialize u at exits from the boundary terminal costs.
       Then BFS outward, computing u[predecessor, exit] from
       u[exit, exit] via the equation
         u[a,b] - u[b,a] = -sign(Δj)·IntegratedMass(Δj, {a,b})
       Validate that the computed u is gauge-consistent
       (i.e. closing any cycle gives ~0 residual). *)
    ...
  ]
```

Then in `BuildMonotoneComparisonData` (or wherever the solution is finalized), call this helper and merge the resulting `u` rules into `solution` *before* passing to `BuildSolverComparisonData`.

### 4.3  Cycle consistency check

The reconstruction is only well-defined if the flow `j` is itself a *valid* MFG equilibrium. If the network has cycles, integrating `u` around any cycle should give 0 (up to numerical tolerance). Add a sanity check that warns if this is violated:

```mathematica
If[Max[cycleResiduals] > 10^-6,
   Message[MonotoneSolver::cyclemismatch, Max[cycleResiduals]];
];
```

### 4.4  Test

Extend `nonlinear-residual.mt` (Phase 1.5) with the same test cases for `MonotoneSolver`. After Phase 4 these should pass.

**Phase 4 exit criteria**:
- `nonlinear-residual.mt` Monotone tests pass.
- `Scripts/VerifyNonLinearSolutions.wls` shows `Mono_NonLinear < 10^-4` for every previously-uncheckable case.
- `MonotoneSolver`'s `AssoMonotone` now contains `u[*]` keys.
- The new cycle-consistency warning fires only on genuinely broken inputs (manual test: feed it a non-Kirchhoff `j` and confirm warning).

---

## Phase 5 — Comprehensive verification (~1 h)

**Goal**: prove the fixes don't regress existing functionality and prove they actually solve the originally-reported problem.

**Actions**:
1. Run `Scripts/RunTests.wls fast` end to end. Must pass.
2. Run `Scripts/RunTests.wls slow` end to end. Must pass (or any new failures are explicitly investigated).
3. Run `Scripts/BenchmarkSuite.wls small` and `Scripts/BenchmarkSuite.wls core`. Must pass with no regressions vs the Phase 0 baseline.
4. Run `Scripts/NonLinearBenchmark.wls`. Every case must show `NonLinearResidual < 10^-4`.
5. Run `Scripts/VerifyNonLinearSolutions.wls`. Every cell in the `NL✓` and `Mono✓` columns must be `YES`.
6. Run `Scripts/CompareDNF.wls --tag "post nonlinear residual fix"` to confirm no DNFReduce regression.

**Phase 5 exit criteria**: all of the above pass; the diff in `DNF_PERFORMANCE_HISTORY.md` shows ≤ 5% slowdown (or is justified by the new correctness check overhead).

**Risk — performance**: computing `NonLinearResidual` per iteration in Phase 3.1 may be expensive on large cases (it `Replace`s into a long symbolic expression). Mitigation: cache the substitution result, or compute only every Nth iteration, or only at convergence-check time. If tests reveal noticeable slowdown, profile and optimize before commit.

---

## Phase 6 — Documentation and changelog (~30 min)

**Goal**: make sure future contributors understand what was wrong and what changed.

**Files to update**:
1. **`CLAUDE.md`** — under "Solver chain & return formats", document the new `NonLinearResidual` field and explain that `KirchhoffResidual` alone is *insufficient* to verify a solution. Add a note in "Debugging & Profiling Workflows" pointing to `Scripts/VerifyNonLinearSolutions.wls`.

2. **`MFGraphs/NonLinearSolver.wl`** — update the `NonLinearSolver::usage` docstring to mention that convergence is now measured against the equation residual, not just `flowDelta`.

3. **`MFGraphs/Monotone.wl`** — update `MonotoneSolver::usage` to mention that `AssoMonotone` now contains reconstructed `u` variables.

4. **New file `BUGFIX_NONLINEAR_RESIDUAL.md`** at repo root — a postmortem-style document explaining:
   - What was broken (with the residual table from Phase 0).
   - Why it went undetected for so long (no test of α≠1).
   - The two-bug structure (preprocessing + convergence check).
   - The fix in each phase.
   - How to verify the fix (`Scripts/VerifyNonLinearSolutions.wls`).
   - How to prevent recurrence (the new regression test in `nonlinear-residual.mt`).

5. **`DNF_PERFORMANCE_HISTORY.md`** and **`PARALLEL_PERFORMANCE_HISTORY.md`** — add an entry tagged `"non-linear residual fix"` with the before/after timings from Phase 5.

6. **Commit messages**: each phase's commit should be self-contained and bisectable. Suggested message format:

   ```
   fix(solver): <one-line summary>

   Phase N of the non-linear residual bug fix (see BUGFIX_NONLINEAR_RESIDUAL.md).

   <one-paragraph explanation of what this commit does and why>

   Verified by:
   - <test or script name>: <result>
   ```

**Phase 6 exit criteria**: documentation is complete, all commits are pushed to the branch, and a draft PR is opened with a checklist matching the phase exit criteria.

---

## Critical files to touch (summary)

| File | Phase | Change |
|------|-------|--------|
| `MFGraphs/MFGraphs.wl` | 1.1, 1.3 | Add `NonLinearResidual` to `BuildSolverComparisonData`; update `IsFeasible`. |
| `MFGraphs/NonLinearSolver.wl` | 1.4, 3.1, 6 | Real `IsNonLinearSolution`; residual-aware convergence; usage docs. |
| `MFGraphs/DataToEquations.wl` | 3.2, 3.3 | Move `u`-equality shortcut out of `MFGPreprocessing`; possibly add gauge anchoring. |
| `MFGraphs/Monotone.wl` | 4.2, 4.3, 6 | Reconstruct `u` post-hoc; cycle consistency check; usage docs. |
| `MFGraphs/Tests/nonlinear-residual.mt` | 1.5 (NEW) | Regression test that fails on master, passes on fix. |
| `Scripts/NonLinearBenchmark.wls` | 1.6 | Surface `NonLinearResidual` in CSV/MD reports. |
| `Scripts/VerifyNonLinearSolutions.wls` | already exists | Use as the manual verification tool throughout. |
| `Scripts/DiagnoseStep*.wls` | 2 (NEW, throwaway) | Phase-2 diagnostic experiments. |
| `CLAUDE.md`, `BUGFIX_NONLINEAR_RESIDUAL.md` | 6 | Documentation. |

## Critical files **not** to touch

- `MFGraphs/DNFReduce.wl` — the DNF engine is innocent here; touching it would risk performance regressions and is unrelated to the bug.
- `MFGraphs/TimeDependentSolver.wl` — out of scope for this fix; the time-dependent solver may share the same bug, but adding it to this PR would balloon the change set. File a follow-up issue.
- `Scripts/BenchmarkSuite.wls` and the benchmark tier definitions — out of scope. Phase 5 *runs* the benchmarks but doesn't modify them.
- The DNF performance history files (until Phase 6 documentation step).

---

## Phase ordering rationale

The order is **deliberate and non-obvious**:

1. **Detection before fix** (Phase 1 before Phase 3). If we fix the convergence check first without adding the residual field, we have no way to verify the fix worked. Worse, the test we add later might pass for the wrong reason.

2. **Diagnosis before fix** (Phase 2 before Phase 3). The two bug hypotheses (preprocessing and convergence) might both be necessary, only one might be necessary, or there might be a third unknown. Fixing prematurely risks chasing the wrong cause and making the code more confusing.

3. **NonLinear before Monotone** (Phase 3 before Phase 4). `NonLinearSolver` is the simpler case (it already has `u` in its solution dict; only the *values* are wrong). `MonotoneSolver` requires structural changes to expose `u`. Fixing the simpler case first builds confidence in the diagnostic infrastructure before tackling the harder one.

4. **Verify before document** (Phase 5 before Phase 6). Documentation that describes a fix that doesn't actually work is worse than no documentation. Phase 5 is the gate.

## Rollback plan

Each phase commits to its own branch tip. If any phase's exit criteria fail and the cause is unclear:

- **Phase 1 rollback**: revert the `IsFeasible` change (the most user-visible). The new fields can stay since they're additive.
- **Phase 2 rollback**: throw away the diagnostic scripts; no production code changed.
- **Phase 3 rollback**: revert the `NonLinearSolver` and `MFGPreprocessing` commits separately. The convergence-check change and the preprocessing change are in different commits *for this exact reason*.
- **Phase 4 rollback**: revert the `Monotone.wl` commit; `nonlinear-residual.mt` Monotone tests stay XFAIL until a follow-up PR.
- **Phase 5/6 rollback**: trivial, doc-only.

## Estimated total effort

| Phase | Time | Risk |
|-------|------|------|
| 0 — Snapshot & branch | 10 min | None |
| 1 — Detection | 2 h | Low (additive) |
| 2 — Diagnose | 1 h | Low (read-only) |
| 3 — Fix NonLinear | 3 h | Medium (logic change) |
| 4 — Fix Monotone | 2 h | Medium (new code path) |
| 5 — Verify | 1 h | Low (running scripts) |
| 6 — Document | 30 min | None |
| **Total** | **~10 h** | |

Plus 1–2 h slack for unexpected discoveries during Phase 2 / 3.

## Open questions for the user before starting

1. **Scope of α=1 fix**: if Phase 2 Diag-4 confirms the linear baseline is *also* affected, do we widen this PR to fix that too, or split it into a follow-up? My recommendation: widen, because the fix is the same code change and a half-fixed solver is worse than the current state.

2. **Approach for Bug A**: Option A (flag) or Option B (move shortcut to `DirectCriticalSolver` only)? My recommendation: Option B — it eliminates the bug at the source, makes `MFGPreprocessing` simpler, and is no harder to implement.

3. **MonotoneSolver scope**: Phase 4 (post-hoc `u` reconstruction) is the right fix, but if you'd rather ship the NonLinear fix first and tackle Monotone in a separate PR, that is also reasonable. The Phase 1 detection infrastructure makes Monotone's status visible (`Missing["MissingVariables"]`) so users at least know not to trust it.

4. **Test execution time**: `nonlinear-residual.mt` will add ~2–5 minutes to the `fast` suite (currently ~27 min). Acceptable, or should the new tests go in `slow`?

5. **PR strategy**: one big PR for all six phases, or one PR per phase? My recommendation: one PR for Phases 0–3 (the core fix), a second PR for Phase 4 (Monotone), a third doc-only PR for Phase 6. Phases 1, 2, 5 are part of the core PR.
