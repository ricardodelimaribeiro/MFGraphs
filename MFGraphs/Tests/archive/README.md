# Archived tests

This directory contains dormant tests that are intentionally excluded from the active suites.

- `kkt-gate.mt` depends on `ClassifyKKT`, which was removed from the public API during the critical-only surface refactor; reactivate only if an equivalent public KKT classification contract returns.
- `numeric-state.mt` exercises internal-only Phase 5/6 Fictitious Play backend symbols; reactivate when those backend components are promoted to supported public API surface (see issues #72-#75).
- `criticalcongestion.mt` is a legacy `TestSuite[...]` aggregator superseded by explicit suite wiring in `Scripts/RunTests.wls`.
- `solve-mfg-legacy.mt` contains legacy `SolveMFG` routing coverage for removed NonLinear/Monotone/NormalizeEdgeFunction public surface.
- `solver-contracts-legacy.mt` contains return-shape contracts for removed NonLinear/Monotone public solvers.

Additional tests archived on April 23, 2026:
- `J9F-critical_congestion.mt`, `Jm9-critical_congestion.mt`, `critical-basic-cases.mt`,
  `critical-numeric-backend.mt`, `exact-mode.mt`, `exit-pass-through.mt`,
  `exit-pass-through-nonzero.mt`, `graphics-public-api.mt`,
  `inconsistent-switching-critical.mt`, `infeasible-status.mt`, `run-tests-smoke.mt`,
  `solve-mfg.mt`, `solver-contracts.mt`, `symbolic-underdetermined.mt`,
  `package-loading.mt`.
- Reason: these suites are not focused on `scenarioTools.wl`, `systemTools.wl`, or `unknownsTools.wl`.
