# MFGraphs Test Layout

This directory mixes active CI tests, long-running paper-scale tests, and legacy coverage retained for reference.

## Runner suites (`Scripts/RunTests.wls`)

- `fast`: active scenario/unknowns/system kernel checks.
- `slow`: active non-kernel checks (currently `graphics.mt`).
- `all`: `fast + slow`.
- `legacy`: `legacy-fast + legacy-slow` (currently empty placeholders).
- `full`: `all + legacy`.

## File groups

- Active root tests:
  - Fast suite: `scenario-kernel.mt`, `make-unknowns.mt`, `reduce-system.mt`, `scenario-consistency.mt`
  - Slow suite: `graphics.mt`
- Archived tests: `archive/` (dormant tests with removed or unsupported surfaces and APIs).

## Notes

- `reference_solutions.json` is generated fixture data for solver regression checks.
- Keep `archive/` tests out of active suites unless their dependent public API is restored.
