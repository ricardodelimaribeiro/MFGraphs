# Scripts Layout

This directory contains active maintenance scripts and an `archive/` folder with legacy or experiment-only scripts.

## Test and docs infrastructure

- `RunTests.wls`: suite runner (`fast`, `all`, `oracle`, `archive`, `full`).
- `RunSingleTest.wls`: run one `.mt` file and print pass/fail summary.
- `GenerateDocs.wls`: regenerate `API_REFERENCE.md` from `::usage` strings; refuses to publish when the usage-arity lint fails.
- `UsageArityLint.wls`: shared usage-arity lint sourced by `GenerateDocs.wls` and `MFGraphs/Tests/usage-arity.mt` (not a standalone CLI).
- `AuditCuratedResults.wls`: recheck exact soundness of the selected saved WXF artifacts and export readable rules, domains, and original constraints into a new directory; does not solve again or replace saved completeness outcomes.
- `DiagnoseArtifactRoundtrip.wls`: read-only diagnosis of original vs WXF-reloaded mathematical content; writes a new evidence directory and never overwrites the failing run.
- `AuditPublicAPI.py`: API surface inventory/audit helper.
- `CheckCriticalSurfaceTests.wls`: gate checks for critical-surface test policy.
- `CheckReportPlaceholders.wls`: fails when committed history/report files still contain placeholder prose.

## Benchmarks

- `BenchmarkSystemSolver.wls`: **the active solver benchmark** (required for before/after runs on solver-sensitive PRs). Flags: `--solver dnf|optimizeddnf|activeset|reduce|boolean|findinstance`, `--tag`, `--timeout`, `--case`. Appends a history entry to `BENCHMARKS.md` and stores solutions under `Results/`.
- `BenchmarkReduceSystem.wls`: historical focused benchmark for raw `reduceSystem`.
- `BenchmarkDNF.wls`: DNF reducer micro-benchmark (assumes cwd is the repo root).
- `BenchmarkBFSDNFReduce.wls`: head-to-head `dnfReduceSystem` (DFS) vs `bfsDNFReduceSystem` (BFS).
- `BenchmarkActiveSetLPPrecheck.wls`: head-to-head `activeSetReduceSystem` with/without `LPPrecheck`.
- `BaselinePruningBenchmark.wls`: measures the paper §5.4 `booleanReduceSystem` baseline (materializes all complementarity branches).
- `PaperBenchmark.wls`: paper §5.4 benchmark (A.2) + C.9 smoke test; modes `merge|fork|jamarat|all`.
- `RepeatJamarat.wls`: repeated timing of the Jamarat critical-congestion cell (the paper's "under five minutes" claim); reports median/mean/min-max over N runs. Flags: `--runs`, `--timeout`, `--memcap`.
- `SweepBranchStateThreshold.wls`: sweeps `$optimizedDNFVarThreshold` to find the branch-state/dnfReduce crossover.

## Profilers

- `ProfileScenarioKernel.wls`: non-mutating staged profiler for scenario → unknowns → system → solver (`--solver`, `--case`, `--timeout`).
- `ProfileDNFReduce.wls`: DNF reducer profiler (`--case name|all`, `--order`, `--timeout`).
- `ProfileDNFProcedural.wls`: baseline profiler for `dnfReduceProcedural`.

## Solver comparison and diagnostics

- `CompareDNFRecursiveProcedural.wls`: correctness comparison of `dnfReduce` (recursive) vs `dnfReduceProcedural` (stack-based); writes a markdown report.
- `CompareSingleCase.wls`: runs the solver columns on one case in a fresh kernel (for big grids); delegates single stages to `CompareSingleCaseBC.wls` (BooleanConvert), `CompareSingleCaseBM.wls` (`booleanMinimizeSystem`), and `CompareSingleCaseBMR.wls` (`booleanMinimizeReduceSystem`).
- `CompareSolverSolutions.wls`: compares residuals/rules between stored benchmark solution files (`Results/*_latest.wl`).
- `AnalyzeSolutionBranches.wls`: analyzes residual solution branches by entry costs.
- `NonCriticalDiagnostic.wls`: convergence diagnostic for the iterative non-critical (Alpha != 1) solver.

## Solution cache

- `RegenerateSolutions.wls`: regenerates cached example solutions under `<repo>/solutions/*.wxf`.
- `SolutionCacheHelpers.wls`: shared cache helpers sourced by `RegenerateSolutions.wls` and `MFGraphs/Tests/example-coverage.mt` (not a standalone CLI).

## Paper-support utilities

`PaperBenchmark.wls` and `BaselinePruningBenchmark.wls`, described above, are
the package-level computation harnesses supporting the stationary-MFG paper.
Manuscript-specific PDF generators live with the manuscript in its `scripts/`
directory, beside the `figures/` and `executables/` folders they produce.
See [`../REPRODUCIBILITY.md`](../REPRODUCIBILITY.md) for the clean-clone
validation and release workflow.

- `ExtractPaperResults.wls`: extracts/validates paper-tier benchmark results (consumes the archived `BenchmarkSuite.wls` output format).
- `MergeBenchmarkResults.py`: merges `Results/benchmark_*.json` outputs (archived pipeline format).

## Research

- `Research/`: self-contained investigation scripts with their committed result data (branching studies and related experiments).

## One-off snapshots

- `perf_review_targeted.wls` + `perf_review_targeted_results.json`: targeted performance-review snapshot.
- `benchmark_preprocessing_results.json`: preprocessing benchmark snapshot (its producer now lives in `archive/`).

## Legacy notes

- `LEGACY_GETEXAMPLEDATA.md`
- `LEGACY_SOLVER_SCRIPTS.md`

## Archive

- `Scripts/archive/` holds older benchmark/profiling/CI scripts not used in the active scenario-kernel workflow.
- `Scripts/repro/` holds preserved debugging and reproduction scripts from active investigation sessions.

## Curated exact research collection

For the bounded paper-foundation collection, use
`python3 Scripts/run_curated_exact.py --help` and the
[technical report](../docs/research/exact-foundation.md).
`CuratedExactScenarios.wls` contains exact inputs. `RunCuratedExact.wls` is the
fresh-kernel worker; `ExactArtifactIntegrity.wls` checks reopened content separately
from `exactSolutionReport`'s mathematical proofs. Every invocation saves a unique
directory under `Results/exact-foundation/`, including timeouts and skipped retries.
Run timing experiments sequentially. For an interleaved source comparison:

```bash
python3 Scripts/run_curated_exact.py --cases grid3,grid4 --methods dnf,linear-net \
  --compare-root PATH_TO_EXTRACTED_BASELINE --baseline-methods dnf \
  --timeout 60 --validation-timeout 5 --repetitions 2 --tag comparison
```

Inspect a saved exact solution and parameter domain without solving again:

```bash
wolframscript -file Scripts/InspectCuratedResult.wls PATH_TO_RESULT_WXF
```
