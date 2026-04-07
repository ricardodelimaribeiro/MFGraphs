# Baseline: Non-linear Residual Bug — Broken State

**Date**: 2026-04-07
**Branch (snapshot)**: `claude/vibrant-aryabhata` → `fix/nonlinear-residual-bug`
**Mathematica**: 14.0.0 for Mac OS X ARM (64-bit)

## Purpose

This document captures the **broken state** of `NonLinearSolver` and `MonotoneSolver`
before any fix is applied. It exists so that:

1. The bug is recorded in git history (not just in conversations).
2. Phase 1 of the fix can produce a regression test that **fails on this commit**
   and **passes after the fix** — proving the fix did something.
3. Future contributors can `git blame` this file to find the postmortem.

See also:
- `.claude/plans/fix-nonlinear-residual-bug.md` — the full multi-phase fix plan
- `Scripts/VerifyNonLinearSolutions.wls` — the script that produced these numbers
- `Scripts/DebugNonLinearResidual.wls` — single-case deep-dive that pinpointed the bug

## What's broken

Both `NonLinearSolver` and `MonotoneSolver` report `KirchhoffResidual = 0` and
`Status = "Feasible"` on every test case below, but their solutions **do not satisfy**
the actual non-linear MFG equation `Nlhs == Nrhs`. The reported `KirchhoffResidual`
only checks linear flow balance `KM·j = B`; the actual equation residual
`|Nlhs - Nrhs|∞` is between 25 and 50 across all cases.

The bug is silent because:
1. `BuildSolverComparisonData` (`MFGraphs/MFGraphs.wl:117-122`) does not compute
   the non-linear residual at all.
2. `IsNonLinearSolution` (`MFGraphs/NonLinearSolver.wl:228-269`) prints the residual
   diagnostically but does not gate on it.
3. The `solver-contracts.mt` test suite only exercises α=1 cases, where the bug
   happens to be invisible.

## Verification table

Produced by running `wolframscript -file Scripts/VerifyNonLinearSolutions.wls`
on this commit:

| Case            | α    | NL Kirchhoff | NL `\|Nlhs−Nrhs\|∞` | Mono Kirchhoff | Mono `\|Nlhs−Nrhs\|∞` | NL OK | Mono OK |
|-----------------|------|--------------|----------------------|-----------------|-------------------------|-------|---------|
| Y7_alpha2       | 2    | 0.           | **50.0**             | 0.              | Missing[NoNumeric]      | NO    | NO      |
| Y7_alpha125     | 1.25 | 0.           | **50.0**             | 0.              | Missing[NoNumeric]      | NO    | NO      |
| Y7_alpha15      | 1.5  | 0.           | **50.0**             | 0.              | Missing[NoNumeric]      | NO    | NO      |
| Case12_alpha25  | 2.5  | 0.           | **25.0**             | 0.              | Missing[NoNumeric]      | NO    | NO      |
| Braess_alpha15  | 1.5  | Missing      | Missing              | Missing         | Missing                 | NO    | NO      |
| Case12_alpha3   | 3    | 0.           | **25.0**             | 0.              | Missing[NoNumeric]      | NO    | NO      |

**Legend**:
- `Kirchhoff` = `Norm[KM·j - B, Infinity]` — linear flow balance, **always 0** (this is what current code reports as the residual)
- `|Nlhs−Nrhs|∞` = `Max[Abs[Nlhs - Nrhs]]` evaluated against the solver's solution — **the actual MFG equation residual**
- `NL/Mono OK` = `True` iff the actual residual is below `10⁻⁴`

## Single-case deep-dive (Y7, α=2)

From `Scripts/DebugNonLinearResidual.wls`:

- **Topology**: Y-network with entrance at vertex 1 (`I1=50`), exits at vertices 3 and 4 (`U1=U2=0`).
- **Solver**: `NonLinearSolver[d2e, "CongestionExponentFunction" -> Function[edge, 2]]`
- **Returned `j` values** (correct):
  ```
  j[en1, 1]    -> 50    (entry boundary, correct)
  j[1, 2]      -> 50    (single path through bottleneck)
  j[2, 3]      -> 25    (split equally)
  j[2, 4]      -> 25    (split equally)
  j[3, ex3]    -> 25
  j[4, ex4]    -> 25
  ```
- **Returned `u` values** (ALL ZERO — the bug):
  ```
  u[1,2] -> 0,  u[2,1] -> 0
  u[2,3] -> 0,  u[3,2] -> 0
  u[2,4] -> 0,  u[4,2] -> 0
  ... (every u variable)
  ```
- **Equation `Nlhs == Nrhs` evaluated against solution**:
  ```
  Nlhs = {50, 25, 25}    (j[a,b]-j[b,a]+u[a,b]-u[b,a])
  Nrhs = {≈0, ≈0, ≈0}    (j[a,b]-j[b,a] - sign·IntegratedMass)
  Diff = {50, 25, 25}    (the equation is grossly violated)
  ```
- The equation reduces to `u[a,b] - u[b,a] = -sign(Δj)·IntegratedMass(Δj, edge)`,
  which for `j[1,2]=50` and `IntegratedMass[50, {1,2}] ≈ 50` requires
  `u[1,2] - u[2,1] = -50`. The solver returns `0`.

## Root cause hypotheses (to be confirmed in Phase 2)

1. **Bug A** — `MFGraphs/DataToEquations.wl:455-470`: `MFGPreprocessing` collapses
   `IneqSwitchingByVertex` to `u`-equalities when all switching costs are zero, then
   `TripleClean` propagates `u = 0` through the entire graph from the boundary
   conditions before `MFGSystemSolver` ever runs.
2. **Bug B** — `MFGraphs/NonLinearSolver.wl:130-136`: the `NestWhileList` convergence
   predicate only checks `Δj`, never the equation residual. The iteration exits as
   `"ToleranceMet"` after one step because the seed already has correct `j`.
3. **Bug C** — `MFGraphs/Monotone.wl`: `AssoMonotone` does not contain `u[*]` keys,
   so the residual is uncheckable.
4. **Bug D** — `MFGraphs/MFGraphs.wl:117-122`: `BuildSolverComparisonData` reports
   only `KirchhoffResidual`, never `NonLinearResidual`.
5. **Bug E** — `MFGraphs/NonLinearSolver.wl:228-269`: `IsNonLinearSolution` is a
   verbose printer that returns the assoc unconditionally, not a true predicate.

## Test coverage gap

`MFGraphs/Tests/solver-contracts.mt` calls `NonLinearSolver` with
`"CongestionExponentFunction" -> Function[edge, 1]` only — i.e., the linear case
where the bug does not visibly trip. **There is no automated test that exercises
α ≠ 1.** This is why the bug has survived for the entire history of the codebase.

The fix plan addresses this in Phase 1.5 with a new file
`MFGraphs/Tests/nonlinear-residual.mt` containing 8 test cases organized along
three orthogonal axes (topology, α regime, switching-cost regime) — see the plan
file for the full matrix.

## Reproduction instructions

```bash
# From repo root, on this commit:
wolframscript -file Scripts/VerifyNonLinearSolutions.wls
# Expect: "NL✓ NO" and "Mono✓ NO" on every row.

# For a single-case deep dive:
wolframscript -file Scripts/DebugNonLinearResidual.wls
# Expect: "Max abs residual: 50.0" and all u variables → 0.
```
