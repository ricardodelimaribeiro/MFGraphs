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
_TODO_

#### Interpretation
_TODO_

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
_TODO_

#### Interpretation
_TODO_

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