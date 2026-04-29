# MFGraphs Test Layout

This directory contains the active scenario-kernel tests and archived legacy/solver-oriented suites.

## Runner suites (`Scripts/RunTests.wls`)

- `fast`: active scenario-kernel suite.
- `all`: alias for `fast`.
- `archive`: archived compatibility/legacy suites (explicit use only).
- `full`: `fast + archive`.

## File groups

- Active tests: root `.mt` files (`scenario-kernel`, `symbolic-unknowns`, `reduce-system`, `scenario-consistency`, `graphicsTools`).
- Archived tests: `archive/` (legacy or non-core surfaces).

## Notes

- `reference_solutions.json` is generated fixture data for solver regression checks.
- Keep `archive/` tests out of CI unless intentionally validating legacy compatibility.
