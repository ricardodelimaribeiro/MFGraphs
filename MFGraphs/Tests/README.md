# MFGraphs Test Layout

This directory mixes active CI tests, long-running paper-scale tests, and legacy coverage retained for reference.

## Runner suites (`Scripts/RunTests.wls`)

- `fast`: critical-only public CI gate.
- `slow`: long-running paper-scale cases.
- `all`: `fast + slow`.
- `meta`: runner-level smoke checks.
- `full`: `all + meta`.
- `unsupported`: legacy/non-critical tests excluded from active CI.

## File groups

- Core/active tests: most `.mt` files at this directory root.
- Meta tests: `run-tests-smoke.mt`.
- Unsupported tests:
  - `inconsistent-switching-critical.mt`
  - `symbolic-underdetermined.mt`
- Archived tests: `archive/` (dormant tests with removed or unsupported surfaces).

## Notes

- `reference_solutions.json` is generated fixture data for solver regression checks.
- Keep `archive/` tests out of active suites unless their dependent public API is restored.
