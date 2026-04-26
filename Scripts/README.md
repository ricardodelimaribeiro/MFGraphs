# Scripts Layout

This directory contains active maintenance scripts and an `archive/` folder with legacy or experiment-only scripts.

## Active scripts

- `RunTests.wls`: suite runner (`fast`, `all`, `archive`, `full`).
- `RunSingleTest.wls`: run one `.mt` file and print pass/fail summary.
- `GenerateDocs.wls`: regenerate API documentation artifacts.
- `AuditPublicAPI.py`: API surface inventory/audit helper.
- `CheckCriticalSurfaceTests.wls`: gate checks for critical-surface test policy.
- `BenchmarkReduceSystem.wls`: focused benchmark for `ReduceSystem`.

## Analysis/report helpers

- `ExtractPaperResults.wls`, `RegenPaperFigures.wls`, `CheckReportPlaceholders.wls`
- `MergeBenchmarkResults.py`

## Legacy notes

- `LEGACY_GETEXAMPLEDATA.md`
- `LEGACY_SOLVER_SCRIPTS.md`

## Archive

- `Scripts/archive/` holds older benchmark/profiling/CI scripts not used in the active scenario-kernel workflow.
