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
| 1 |           -6
3.000 × 10 |           -6
8.000 × 10 |           -6
2.000 × 10 | 0.00001 | OK |
| 2 |           -6
3.000 × 10 |           -6
5.000 × 10 |           -6
4.000 × 10 |           -6
9.000 × 10 | OK |
| 3 |           -6
2.000 × 10 |           -6
3.000 × 10 |           -6
6.000 × 10 |           -6
9.000 × 10 | OK |
| 4 |           -6
3.000 × 10 |           -6
3.000 × 10 |           -6
4.000 × 10 |           -6
7.000 × 10 | OK |
| 5 |           -6
4.000 × 10 |           -6
3.000 × 10 |           -6
3.000 × 10 |           -6
6.000 × 10 | OK |
| 6 |           -6
3.000 × 10 |           -6
4.000 × 10 |           -6
9.000 × 10 | 0.000 | OK |
| 27 |           -6
5.000 × 10 |           -6
3.000 × 10 |           -6
7.000 × 10 | 0.00001 | OK |

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
| 1 |           -6
4.000 × 10 |           -6
5.000 × 10 |           -6
3.000 × 10 |           -6
8.000 × 10 | OK |
| 2 |           -6
5.000 × 10 |           -6
3.000 × 10 |           -6
9.000 × 10 | 0.000 | OK |
| 3 |           -6
3.000 × 10 |           -6
3.000 × 10 | 0.000 | 0.000 | OK |
| 4 |           -6
4.000 × 10 |           -6
6.000 × 10 |           -6
6.000 × 10 | 0.000 | OK |
| 5 |           -6
4.000 × 10 |           -6
4.000 × 10 |           -6
7.000 × 10 | 0.000 | OK |
| 6 |           -6
3.000 × 10 |           -6
3.000 × 10 | 0.000 | 0.000 | OK |
| 27 |           -6
7.000 × 10 |           -6
4.000 × 10 | 0.000 | 0.000 | OK |

#### Rationale
_baseline serial_

#### Changes
_TODO_

#### Interpretation
_TODO_

---

### 2026-03-29 — baseline serial

**Commit:** `188e002`
**Date:** Sun 29 Mar 2026 15:52:12
**Kernels:** 16

#### Tier: small

| Case | D2E (s) | CritSolver (s) | NLSolver (s) | Total (s) | Status |
|------|---------|---------------|-------------|-----------|--------|
| 1 | 0. | 0. | 0. | 0. | OK |
| 2 | 0. | 0. | 0. | 0. | OK |
| 3 | 0. | 0. | 0. | 0. | OK |
| 4 | 0. | 0. | 0. | 0. | OK |
| 5 | 0. | 0. | 0. | 0. | OK |
| 6 | 0. | 0. | 0. | 0. | OK |
| 27 | 0. | 0. | 0. | 0. | OK |

#### Rationale
_baseline serial_

#### Changes
_TODO_

#### Interpretation
_TODO_

---

### 2026-03-31 — Phase A+B+C parallel

**Commit:** `0015297`
**Date:** Tue 31 Mar 2026 12:05:01
**Kernels:** 16

#### Tier: small

| Case | D2E (s) | CritSolver (s) | NLSolver (s) | Total (s) | Status |
|------|---------|---------------|-------------|-----------|--------|
| 1 | 0. | 0. | 0. | 0. | OK |
| 2 | 0. | 0. | 0. | 0. | OK |
| 3 | 0. | 0. | 0. | 0. | OK |
| 4 | 0. | 0. | 0. | 0. | OK |
| 5 | 0. | 0. | 0. | 0. | OK |
| 6 | 0. | 0. | 0. | 0. | OK |
| 27 | 0. | 0. | 0. | 0. | OK |

#### Rationale
_Phase A+B+C parallel_

#### Changes
_TODO_

#### Interpretation
_TODO_

---

### 2026-03-31 — post-merge master

**Commit:** `6fd0415`
**Date:** Tue 31 Mar 2026 12:10:16
**Kernels:** 16

#### Tier: small

| Case | D2E (s) | CritSolver (s) | NLSolver (s) | Total (s) | Status |
|------|---------|---------------|-------------|-----------|--------|
| 1 | 0. | 0. | 0. | 0. | OK |
| 2 | 0. | 0. | 0. | 0. | OK |
| 3 | 0. | 0. | 0. | 0. | OK |
| 4 | 0. | 0. | 0. | 0. | OK |
| 5 | 0. | 0. | 0. | 0. | OK |
| 6 | 0. | 0. | 0. | 0. | OK |
| 27 | 0. | 0. | 0. | 0. | OK |

#### Rationale
_post-merge master_

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