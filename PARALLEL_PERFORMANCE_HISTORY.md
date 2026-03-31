# Parallel Solver Performance History

Records wall-clock timing impact of each parallelization change to the MFGraphs solver
pipeline. Entries are usually generated automatically by `BenchmarkSuite.wls --tag`,
with occasional manual notes when a benchmark run stalls or needs interpretation.

**How to add an entry:**
```bash
wolframscript -file Scripts/BenchmarkSuite.wls small tag="description of change"
wolframscript -file Scripts/BenchmarkSuite.wls core tag="description of change"
wolframscript -file Scripts/BenchmarkSuite.wls stress tag="description of change"
wolframscript -file Scripts/BenchmarkSuite.wls stress case=7 timeout=3600
```

**Interpreting results:** Speedup is only meaningful relative to an entry with the same
number of kernels (`$KernelCount`). A baseline entry (`$KernelCount = 0`) records
serial timings against which parallel entries can be compared.

Entries generated before the current benchmark schema may omit the `Monotone (s)`
column and may report a narrower status summary. Newer entries include all three
stationary solver times plus a case-level status that surfaces any solver-specific
timeouts or failures.

---

### 2026-03-29 — baseline serial

**Commit:** `199012c`
**Date:** Sun 29 Mar 2026 15:36:17
**Kernels:** 16

#### Tier: small

| Case | D2E (s) | CritSolver (s) | NLSolver (s) | Total (s) | Status |
|------|---------|---------------|-------------|-----------|--------|
| 1 | 3.000 × 10^-6 | 8.000 × 10^-6 | 2.000 × 10^-6 | 0.00001 | OK |
| 2 | 3.000 × 10^-6 | 5.000 × 10^-6 | 4.000 × 10^-6 | 9.000 × 10^-6 | OK |
| 3 | 2.000 × 10^-6 | 3.000 × 10^-6 | 6.000 × 10^-6 | 9.000 × 10^-6 | OK |
| 4 | 3.000 × 10^-6 | 3.000 × 10^-6 | 4.000 × 10^-6 | 7.000 × 10^-6 | OK |
| 5 | 4.000 × 10^-6 | 3.000 × 10^-6 | 3.000 × 10^-6 | 6.000 × 10^-6 | OK |
| 6 | 3.000 × 10^-6 | 4.000 × 10^-6 | 9.000 × 10^-6 | 0.000 | OK |
| 27 | 5.000 × 10^-6 | 3.000 × 10^-6 | 7.000 × 10^-6 | 0.00001 | OK |

#### Rationale
_baseline serial_

#### Changes
_Initial benchmark checkpoint captured before later solver and reporting changes added the monotone column and case-level status summary._

#### Interpretation
_All recorded small-tier cases completed essentially instantaneously. Treat the exact totals in this legacy entry as historical context only, because they predate the newer schema that records monotone timings explicitly._

---

### 2026-03-29 — baseline serial

**Commit:** `188e002`
**Date:** Sun 29 Mar 2026 15:44:02
**Kernels:** 16

#### Tier: small

| Case | D2E (s) | CritSolver (s) | NLSolver (s) | Total (s) | Status |
|------|---------|---------------|-------------|-----------|--------|
| 1 | 4.000 × 10^-6 | 5.000 × 10^-6 | 3.000 × 10^-6 | 8.000 × 10^-6 | OK |
| 2 | 5.000 × 10^-6 | 3.000 × 10^-6 | 9.000 × 10^-6 | 0.000 | OK |
| 3 | 3.000 × 10^-6 | 3.000 × 10^-6 | 0.000 | 0.000 | OK |
| 4 | 4.000 × 10^-6 | 6.000 × 10^-6 | 6.000 × 10^-6 | 0.000 | OK |
| 5 | 4.000 × 10^-6 | 4.000 × 10^-6 | 7.000 × 10^-6 | 0.000 | OK |
| 6 | 3.000 × 10^-6 | 3.000 × 10^-6 | 0.000 | 0.000 | OK |
| 27 | 7.000 × 10^-6 | 4.000 × 10^-6 | 0.000 | 0.000 | OK |

#### Rationale
_baseline serial_

#### Changes
_Initial baseline recorded before parallelization changes._

#### Interpretation
_All small-tier cases complete in sub-millisecond time._



### 2026-03-31 — post-simplify baseline

**Commit:** `78114db`
**Date:** Tue 31 Mar 2026 15:20:56
**Kernels:** 16

#### Tier: small

| Case | D2E (s) | CritSolver (s) | NLSolver (s) | Total (s) | Status |
|------|---------|---------------|-------------|-----------|--------|
| 1 | 3.e-6 | 4.e-6 | 2.e-6 | 6.e-6 | OK |
| 2 | 4.e-6 | 4.e-6 | 3.e-6 | 7.e-6 | OK |
| 3 | 3.e-6 | 5.e-6 | 4.e-6 | 9.e-6 | OK |
| 4 | 4.e-6 | 3.e-6 | 4.e-6 | 7.e-6 | OK |
| 5 | 4.e-6 | 6.e-6 | 3.e-6 | 9.e-6 | OK |
| 6 | 4.e-6 | 3.e-6 | 5.e-6 | 8.000000000000001e-6 | OK |
| 27 | 6.e-6 | 4.e-6 | 4.e-6 | 8.e-6 | OK |

#### Rationale
_post-simplify baseline_

#### Changes
_Recorded after the simplify follow-up work to confirm that the small-tier baseline stayed stable before the later benchmark schema refresh._

#### Interpretation
_Small-tier timings stayed within the same microsecond-scale band as the earlier baseline, with no visible regression in the critical or nonlinear stationary solvers._

---

### 2026-04-01 — post-merge small-tier sanity baseline

**Commit:** `e59a562`
**Date:** Wed 1 Apr 2026 01:23:17
**Kernels:** 0

#### Tier: small

| Case | D2E (s) | CritSolver (s) | NLSolver (s) | Monotone (s) | Total (s) | Case Status |
|------|---------|---------------|-------------|--------------|-----------|-------------|
| 1 | 4.e-6 | 5.e-6 | 3.e-6 | 1.e-5 | 1.8e-5 | OK |
| 2 | 4.e-6 | 4.e-6 | 5.e-6 | 1.e-5 | 1.9e-5 | OK |
| 3 | 6.e-6 | 4.e-6 | 4.e-6 | 7.e-6 | 1.5e-5 | OK |
| 4 | 4.e-6 | 6.e-6 | 4.e-6 | 4.e-6 | 1.4e-5 | OK |
| 5 | 5.e-6 | 5.e-6 | 5.e-6 | 6.e-6 | 1.6e-5 | OK |
| 6 | 6.e-6 | 8.e-6 | 5.e-6 | 5.e-6 | 1.8e-5 | OK |
| 27 | 6.e-6 | 7.e-6 | 7.e-6 | 6.e-6 | 2.e-5 | OK |

#### Rationale
_Confirm the current `master` branch still has a clean small-tier baseline after the CI and solver changes were merged._

#### Changes
_No new parallelization change in this entry. This is a fresh benchmark checkpoint on merged `master`, recorded alongside a stress-tier investigation._

#### Interpretation
_Small-tier performance remains effectively unchanged and all 21 solver-case runs completed successfully. Attempts to replace the old timeout observations with real numbers show that the slow path is genuinely long: a targeted rerun on the same commit using `wolframscript -file Scripts/BenchmarkSuite.wls stress case=7 timeout=3600` stayed CPU-bound in the `Monotone` path for more than 12 minutes without finishing, and an earlier untargeted medium run remained active for more than 30 minutes. Treat stress-tier comparisons as unreliable until timeout enforcement is fixed or the slow path is isolated._

---

## Future entries (template)

```markdown
### YYYY-MM-DD — Short description

**Commit:** `<hash>` — *"Commit subject"*
**Date:** <date>
**Kernels:** <$KernelCount>

#### Tier: small

| Case | D2E (s) | CritSolver (s) | NLSolver (s) | Monotone (s) | Total (s) | Case Status |
|------|---------|---------------|-------------|--------------|-----------|-------------|
| ...  | ...     | ...           | ...         | ...          | ...       | ...         |

#### Tier: core

| Case | D2E (s) | CritSolver (s) | NLSolver (s) | Monotone (s) | Total (s) | Case Status |
|------|---------|---------------|-------------|--------------|-----------|-------------|
| ...  | ...     | ...           | ...         | ...          | ...       | ...         |

#### Tier: stress

| Case | D2E (s) | CritSolver (s) | NLSolver (s) | Monotone (s) | Total (s) | Case Status |
|------|---------|---------------|-------------|--------------|-----------|-------------|
| ...  | ...     | ...           | ...         | ...          | ...       | ...         |

#### Rationale
_Why was this change made?_

#### Changes
_What was parallelized._

#### Interpretation
_Which cases improved/regressed and why._
```
