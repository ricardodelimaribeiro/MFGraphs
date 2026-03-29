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
