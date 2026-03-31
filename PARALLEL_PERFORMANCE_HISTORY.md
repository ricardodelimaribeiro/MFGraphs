# Parallel Solver Performance History

Records wall-clock timing impact of each parallelization change to the MFGraphs solver
pipeline. Each entry is generated automatically by `BenchmarkSuite.wls --tag`.

**How to add an entry:**
```bash
wolframscript -file Scripts/BenchmarkSuite.wls small  --tag "description of change"
wolframscript -file Scripts/BenchmarkSuite.wls medium --tag "description of change"
```

**Interpreting results:** Speedup is only meaningful relative to an entry with the same
number of kernels (`$KernelCount`). A baseline entry (`$KernelCount = 0`) records
serial timings against which parallel entries can be compared.

---

### 2026-03-29 — baseline serial

**Commit:** `199012c`
**Date:** Sun 29 Mar 2026 15:36:17
**Kernels:** 16

#### Tier: small

| Case | D2E (s) | CritSolver (s) | NLSolver (s) | Total (s) | Status |
|------|---------|---------------|-------------|-----------|--------|
| 1 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | OK |
| 2 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | OK |
| 3 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | OK |
| 4 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | OK |
| 5 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | OK |
| 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | OK |
| 27 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | OK |

> [!NOTE]
> The original entries recorded sub-microsecond times that were corrupted by Mathematica's
> multiline scientific notation formatting. The timing precision was lost; values shown as
> `0.0000` represent times under 0.001s.

#### Rationale
_baseline serial_

#### Changes
_Initial baseline recorded before parallelization changes._

#### Interpretation
_All small-tier cases complete in sub-millisecond time._

---

### 2026-03-31 — Phase A+B+C parallel

**Commit:** `0015297`
**Date:** Tue 31 Mar 2026 12:05:01
**Kernels:** 16

#### Tier: small

| Case | D2E (s) | CritSolver (s) | NLSolver (s) | Total (s) | Status |
|------|---------|---------------|-------------|-----------|--------|
| 1 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | OK |
| 2 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | OK |
| 3 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | OK |
| 4 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | OK |
| 5 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | OK |
| 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | OK |
| 27 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | OK |

#### Rationale
_Phase A+B+C parallel_

#### Changes
_Parallelized IneqSwitching simplification, ineqsByTransition selection, and ineqsByTransition simplification in MFGSystemSolver._

#### Interpretation
_Small-tier cases are too small to benefit from parallelization (below $MFGraphsParallelThreshold). No regression._

---

### 2026-03-31 — post-merge master

**Commit:** `6fd0415`
**Date:** Tue 31 Mar 2026 12:10:16
**Kernels:** 16

#### Tier: small

| Case | D2E (s) | CritSolver (s) | NLSolver (s) | Total (s) | Status |
|------|---------|---------------|-------------|-----------|--------|
| 1 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | OK |
| 2 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | OK |
| 3 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | OK |
| 4 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | OK |
| 5 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | OK |
| 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | OK |
| 27 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | OK |

#### Rationale
_post-merge master_

#### Changes
_Merged parallelize branch into master. All parallel dispatch phases active._

#### Interpretation
_No regression on small-tier cases after merge._

---

## Future entries (template)

```markdown
### YYYY-MM-DD — Short description

**Commit:** `<hash>` — *"Commit subject"*
**Date:** <date>
**Kernels:** <$KernelCount>

#### Tier: small

| Case | D2E (s) | CritSolver (s) | NLSolver (s) | Total (s) | Status |
|------|---------|---------------|-------------|-----------|--------|
| ...  | ...     | ...           | ...         | ...       | ...    |

#### Tier: medium

| Case | D2E (s) | CritSolver (s) | NLSolver (s) | Total (s) | Status |
|------|---------|---------------|-------------|-----------|--------|
| ...  | ...     | ...           | ...         | ...       | ...    |

#### Rationale
_Why was this change made?_

#### Changes
_What was parallelized._

#### Interpretation
_Which cases improved/regressed and why._
```