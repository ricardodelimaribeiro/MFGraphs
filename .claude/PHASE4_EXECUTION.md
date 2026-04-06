# Phase 4: Jamaratv9 Execution & Table Assembly

## Overview
Once Phase 3 validation passes, run the full instrumented Jamaratv9 solver, collect statistics, and generate the final tables for §5.3.

**Expected runtime**: ~12 minutes (same as uninstrumented solver, assuming < 10% overhead verified in Phase 3)

## Step 1: Run the Jamaratv9 solver

From the repo root:

```bash
wolframscript -file Scripts/InstrumentJamarat.wls > Results/jamaratv9_run.log 2>&1 &
```

Or in Mathematica directly:

```mathematica
Get["Scripts/InstrumentJamarat.wls"]
(* This will execute and export prunStats to Results/jamaratv9_pruning_stats.json *)
```

Monitor progress:
```bash
tail -f Results/jamaratv9_run.log
```

## Step 2: Check output

After the run completes, verify:

```bash
# Check that the JSON was exported
ls -lh Results/jamaratv9_pruning_stats.json

# View the statistics
cat Results/jamaratv9_pruning_stats.json | jq '.'
```

Expected JSON structure:
```json
{
  "timestamp": "...",
  "network": "Jamaratv9",
  "q": <measured value>,
  "branchesEntered": <count>,
  "branchesPrunedImmediate": <count>,
  "branchesPrunedElimination": <count>,
  "leavesSurviving": <count>,
  "totalWallTime": <seconds>,
  "passTimings": [
    {"elapsed": <time>, "timestamp": "..."},
    ...
  ],
  "dnfReduceStats": {"branchesIn": <n>, "branchesOut": <m>},
  "reduceDisjunctsStats": {"input": <n>, "output": <m>}
}
```

## Step 3: Generate LaTeX tables from statistics

Create `Scripts/GenerateStatsTables.wls`:

```mathematica
#!/usr/bin/env wolframscript
(* Generate LaTeX tables from jamaratv9_pruning_stats.json *)

(* Load JSON *)
statsPath = FileNameJoin[{Directory[], "Results", "jamaratv9_pruning_stats.json"}];
stats = Import[statsPath, "JSON"];

(* Extract key metrics *)
q = stats["q"];
totalTime = stats["totalWallTime"];
branchesEntered = stats["branchesEntered"];
prunedImm = stats["branchesPrunedImmediate"];
prunedElim = stats["branchesPrunedElimination"];
leavesSurviving = stats["leavesSurviving"];
passCount = Length[stats["passTimings"]];
passTimings = stats["passTimings"];

Print[""];
Print["=== STATISTICS SUMMARY ==="];
Print["Network: Jamaratv9"];
Print["q (uncertain triples): ", q];
Print["Total alternatives: 2^", q];
Print[""];
Print["Branches:"];
Print["  Entered: ", branchesEntered];
Print["  Pruned (immediate): ", prunedImm];
Print["  Pruned (elimination): ", prunedElim];
Print["  Surviving: ", leavesSurviving];
Print[""];
Print["Timing:"];
Print["  Total wall time: ", N[totalTime, 5], " seconds"];
Print["  TripleClean passes: ", passCount];
Print["  Avg per pass: ", N[totalTime / passCount, 5], " seconds"];
Print[""];

(* Generate LaTeX table *)
latex = "% Pruning Statistics Table for §5.3\n\n";
latex = latex <> "\\begin{table}[h!]\n";
latex = latex <> "\\caption{Jamaratv9 solver pruning statistics}\n";
latex = latex <> "\\label{tab:jamaratv9-pruning}\n";
latex = latex <> "\\centering\n";
latex = latex <> "\\begin{tabular}{lrr}\n";
latex = latex <> "\\toprule\n";
latex = latex <> "\\textbf{Metric} & \\textbf{Count} & \\textbf{Percentage} \\\\\n";
latex = latex <> "\\midrule\n";

(* Calculate percentages *)
pctImm = If[branchesEntered > 0, 100 * prunedImm / branchesEntered, 0];
pctElim = If[branchesEntered > 0, 100 * prunedElim / branchesEntered, 0];
pctSurv = If[branchesEntered > 0, 100 * leavesSurviving / branchesEntered, 0];

latex = latex <> "Sign alternatives & $2^{" <> ToString[q] <> "}$ & --- \\\\\n";
latex = latex <> "Branches explored & " <> ToString[branchesEntered] <> " & $100\\%$ \\\\\n";
latex = latex <> "Pruned (immediate) & " <> ToString[prunedImm] <> " & $" <> ToString[N[pctImm, 2]] <> "\\%$ \\\\\n";
latex = latex <> "Pruned (elimination) & " <> ToString[prunedElim] <> " & $" <> ToString[N[pctElim, 2]] <> "\\%$ \\\\\n";
latex = latex <> "Surviving leaves & " <> ToString[leavesSurviving] <> " & $" <> ToString[N[pctSurv, 2]] <> "\\%$ \\\\\n";
latex = latex <> "\\midrule\n";
latex = latex <> "Total solver time & --- & " <> ToString[N[totalTime, 2]] <> " s \\\\\n";
latex = latex <> "TripleClean passes & " <> ToString[passCount] <> " & --- \\\\\n";
latex = latex <> "\\bottomrule\n";
latex = latex <> "\\end{tabular}\n";
latex = latex <> "\\end{table}\n";

(* Save to file *)
tablePath = FileNameJoin[{Directory[], "Results", "jamaratv9_stats_table.tex"}];
Export[tablePath, latex, "Text"];
Print["LaTeX table saved to: ", tablePath];
Print[""];
Print["Copy the following table into §5.3:"];
Print[""];
Print[latex];
```

Run it:
```bash
wolframscript -file Scripts/GenerateStatsTables.wls
```

## Step 4: Assemble tables into main.tex

### Deliverable 1: Parameter table (from Phase 1)

The parameter table is already prepared in `.claude/jamaratv9_latex_tables.tex`. Locate §5.3 in `main.tex` (or wherever the Jamaratv9 discussion is) and insert:

```latex
% === JAMARATV9 NETWORK PARAMETERS ===
% Tables 1-4 from .claude/jamaratv9_latex_tables.tex
```

Copy the four tables from that file.

### Deliverable 2: Pruning statistics table

From the LaTeX generated in Step 3, insert into §5.3 after the parameter table:

```latex
% === JAMARATV9 PRUNING STATISTICS ===
% Generated from Results/jamaratv9_pruning_stats.json
```

Paste the LaTeX table.

### Deliverable 3: Update abstract/conclusion

Replace placeholders with measured values:

**In Abstract:**
- "~12 minutes" → actual measured `totalWallTime` (e.g., "621 seconds" or "~10.4 minutes")
- "~$2^{95}$ alternatives" → if $q \neq 95$, update to $2^q$
- "~M surviving leaves" → actual `leavesSurviving` count

**In Conclusion:**
- Same updates as abstract, plus optionally add:
  - "X% immediately pruned" (use pctImm)
  - "Y% eliminated via symbolic solving" (use pctElim)
  - "Z% surviving" (use pctSurv)

## Step 5: Verify consistency

Check that the counts make sense:

```
branchesEntered ≥ branchesPrunedImmediate + branchesPrunedElimination + leavesSurviving
```

Or:
```
branchesEntered - prunedImm - prunedElim = leavesSurviving
```

If this doesn't hold, the instrumentation has a bug (likely double-counting in DNFReduce vs ReduceDisjuncts).

## Step 6: Create final commit

```bash
git add Results/jamaratv9_pruning_stats.json Results/jamaratv9_run.log
git add Scripts/InstrumentJamarat.wls Scripts/GenerateStatsTables.wls
git commit -m "Add Jamaratv9 pruning statistics and instrumentation

Extracted network parameters and instrumented MFGraphs solver to measure
branch pruning on Jamaratv9 network (2^95 alternatives).

Results:
  - q uncertain triples: <value>
  - Total branches explored: <value>
  - Branches surviving: <value>
  - Wall time: <value> seconds

Tables inserted into §5.3 with parameter and statistics summaries."
```

## If Something Goes Wrong

**JSON import fails**: Check that `Results/jamaratv9_pruning_stats.json` is valid:
```bash
jq . Results/jamaratv9_pruning_stats.json
```

**Counts don't match**: Likely DNFReduce and ReduceDisjuncts are both counting the same branches. Either:
1. Remove the DNFReduce wrapper (let only ReduceDisjuncts count), or
2. Adjust the counter logic to avoid double-counting

**Wall time is much longer than expected**: Instrumentation overhead is > 50%. Switch to `Sow`/`Reap`:

```mathematica
(* In InstrumentJamarat.wls, replace AppendTo with *)
{time, r} = AbsoluteTiming[origTripleClean[args]];
Sow[<|"elapsed" -> t, "timestamp" -> Now|>, "pass"];
```

Then collect with:
```mathematica
{result, {passes}} = Reap[(*solver call*), "pass"];
prunStats["passTimings"] = passes;
```
