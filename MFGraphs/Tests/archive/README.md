# Archived tests

This directory contains dormant tests that are intentionally excluded from the active suites.

- `kkt-gate.mt` depends on `ClassifyKKT`, which was removed from the public API during the critical-only surface refactor; reactivate only if an equivalent public KKT classification contract returns.
- `numeric-state.mt` exercises internal-only Phase 5/6 Fictitious Play backend symbols; reactivate when those backend components are promoted to supported public API surface (see issues #72-#75).
