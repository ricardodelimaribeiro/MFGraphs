# MFGraphs Test Layout

This directory contains the active scenario-kernel tests and archived legacy/solver-oriented suites.

## Runner suites (`Scripts/RunTests.wls`)

- `fast`: the 13 active suites (see below).
- `all`: alias for `fast`.
- `oracle`: `numeric-oracle.mt` + `fictitiousPlay.mt` (the opt-in `numericOracle`` classifier suites).
- `archive`: archived compatibility/legacy suites (explicit use only; they target the removed legacy API and are not expected to pass).
- `full`: `fast` + `fictitiousPlay.mt` + `archive`.

## File groups

- Active tests (`fast`): `scenario-kernel`, `symbolic-unknowns`, `reduce-system`, `scenario-consistency`, `graphicsTools`, `orchestration`, `dnf-reducer`, `boolean-minimize`, `tawaf`, `example-coverage`, `numeric-oracle`, `utilities`, `usage-arity`.
- Oracle extra: `fictitiousPlay` (runs in `oracle` and `full`, not `fast`).
- Archived tests: `archive/` (legacy or non-core surfaces).

## Notes

- The live solution cache used by `example-coverage.mt` is `<repo>/solutions/*.wxf` (see `Scripts/SolutionCacheHelpers.wls`); `reference_solutions.json` is a legacy fixture referenced only by archived scripts.
- Keep `archive/` tests out of CI unless intentionally validating legacy compatibility.
