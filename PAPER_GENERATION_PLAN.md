# Paper Generation Plan (Post-Validation)

## Overview

Once the **paper tier benchmark (HRF Scenario 1)** completes validation, the following steps will generate publication-ready tables and figures for the linear critical congestion paper.

---

## Step 1: Extract & Validate Results (Automated)

**Script**: `Scripts/ExtractPaperResults.wls`

```bash
wolframscript -file Scripts/ExtractPaperResults.wls
```

**Output**:
- ✅ HRF Scenario 1 validation status (must show: Status OK, ResultKind Success, Feasibility Feasible)
- Wall-clock time for CriticalCongestion solver
- Memory usage
- Feasibility confirmation

**Gate**: Must pass before proceeding to paper generation.

---

## Step 2: Aggregate Benchmark Results

**Script**: `Scripts/MergeBenchmarkResults.py` (already in repo)

Merge all tier results into unified dataset:
- `small` tier (cases 1–6, 27): baseline performance
- `core` tier (if run): stable reference cases
- `paper` tier (HRF Scenario 1): publication benchmark

**Output**: Consolidated CSV/JSON for analysis.

---

## Step 3: Generate Paper Tables

**Script**: `Scripts/GeneratePaperTables.py` (already in repo)

Produces:
- **Table 1**: Small tier cases — solver time comparison (CriticalCongestion vs NonLinear vs Monotone)
- **Table 2**: HRF Scenario 1 — detailed equilibrium properties and convergence metrics
- **Table 3** (optional): Scalability across tiers (if large/vlarge tiers run successfully)

**Format**: LaTeX tables suitable for direct inclusion in paper.

**Output Files**:
- `Results/paper_table_*.tex` — LaTeX table snippets
- `Results/paper_table_side_by_side.tex` — Comparison views (already in repo)

---

## Step 4: Validate Reference Matches

**Requirement**: Compare CriticalCongestion output against known exact equilibrium for HRF Scenario 1.

**Status**: Reference hashes stored in `GenerateReferenceHashes.wls` output (if previously computed).

**Action**: If reference exists, validate that flow vector matches algebraically (within machine epsilon).

If no reference: Use this run as the canonical reference for future validation.

---

## Step 5: Prepare Paper Narrative

**Key claims to support with benchmark data:**

1. **Correctness**: HRF Scenario 1 equilibrium is exact (algebraic match with symbolic solver path)
2. **Performance**: CriticalCongestion solves HRF in ~122 seconds (acceptable for reproducibility)
3. **Robustness**: All small-tier cases pass without regression
4. **Scalability** (if tested): Large-tier cases solve within timeout budgets

**Tables to include**:
- Wall-clock times by solver
- Memory usage trends
- Convergence metrics (residuals, iterations)
- Feasibility verification

---

## Timeline

| Phase | Status | Action |
|-------|--------|--------|
| Benchmark small tier | ✅ Complete (21/21 OK) | Reference baseline established |
| Benchmark paper tier | 🔄 Running | Wait for Monotone to finish |
| Extract & validate | ⏳ Pending | Run `ExtractPaperResults.wls` |
| Merge results | ⏳ Pending | Run `MergeBenchmarkResults.py` |
| Generate tables | ⏳ Pending | Run `GeneratePaperTables.py` |
| Paper integration | ⏳ Pending | Embed LaTeX tables in document |

---

## Success Criteria

✅ **Paper tier validates**: HRF Scenario 1 shows Feasible, Success result  
✅ **All tables generate**: LaTeX output is compilable and readable  
✅ **Reference matches**: Flow vectors match exact equilibrium (if available)  
✅ **No regressions**: Small tier baseline unchanged from Phase 5/6 changes  

---

## Notes

- **Monotone solver**: May show "NonConverged" on HRF Scenario; this is **expected** and **acceptable** — Monotone is approximate (ODE-based), while CriticalCongestion is exact (symbolic).
- **Paper scope**: Focus on linear critical (α=1) only. Phase 5/6 (Fictitious Play) are backup infrastructure for future non-linear work.
- **Phase 5/6 not in paper**: These are internal-only v1 due to known vulnerabilities (documented in Issues #72–#75). Do not cite or mention in paper submission.

---

**Next Action**: Monitor paper tier completion, then execute Step 1 (ExtractPaperResults.wls).
