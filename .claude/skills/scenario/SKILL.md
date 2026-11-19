---
name: scenario
description: Walk one MFGraphs scenario through the solver pipeline gate by gate (scenario, unknowns, system, solve, valid, exact) with Scripts/WalkScenario.wls. Use to solve or debug a built-in example or a custom network.
---

Run one scenario through the six **gates** of the reference pipeline and report where it stands. `Scripts/WalkScenario.wls` does the walk; you choose the input, run it, and act on the first gate that did not pass.

## Steps

### 1. Fix the input

- **Built-in example:** take the key from the request, or list the keys:
  ```bash
  wolframscript -file Scripts/WalkScenario.wls --list
  ```
- **Custom network:** write one Wolfram expression into a `.wl` file outside the repo. Its value is a scenario, such as `gridScenario[{3}, {{1, 10}}, {{3, 0}}]` or `getExampleScenario["Grid0303", entries, exits]`, or a raw association for `makeScenario`. The constructors (`gridScenario`, `cycleScenario`, `graphScenario`, `amScenario`, `makeScenario`) are specified in `API_REFERENCE.md`.

Done when you hold a key that `--list` prints, or the path of that file.

### 2. Run the walk

```bash
wolframscript -file Scripts/WalkScenario.wls --case <key>
wolframscript -file Scripts/WalkScenario.wls --scenario-file <path.wl>
```

`--solve-timeout <seconds>` (default 60) bounds the solve gate. `--completeness <seconds>` adds the exact completeness search.

Done when the output ends with a `RESULT` line.

### 3. Report

Give the user every gate line as printed, then the `RESULT`. Quote `ResultKind`, `Soundness` and `Completeness` verbatim: `Completeness NotChecked` and `Completeness Timeout` both mean completeness is not established, and the report says exactly that.

### 4. Take the branch of the first gate that did not pass

| Gate and status | Meaning | Next action |
|---|---|---|
| `scenario FAIL` | The input broke a validation rule; the line names the first one. | Fix the input against the `makeScenario` entry in `API_REFERENCE.md`, then rerun. |
| `system FAIL` | The scenario is outside what the solvers take; the line says why (an edge with `Alpha != 1`, or a switching cost the system cannot use). | Stop and report it as the model scope described in `CLAUDE.md`. |
| `solve TIMEOUT` | The reference solver ran out of time. Feasibility is unknown. | Rerun with a larger `--solve-timeout`; if it still times out, profile with `Scripts/ProfileDNFReduce.wls`. |
| `solve NOSOLUTION` | The reference solver reduced the system to `False`. | Report it as the result. |
| `unknowns FAIL`, `solve FAIL`, `valid FAIL`, `exact FAIL` | A package defect, not an input problem. | Stop. Keep the input file and compare against the open `bug` issues before changing any code. |
| output ends before `RESULT` | The kernel ended the evaluation (recursion limit or crash). The exit code is still 0. | Report the walk as incomplete at the last gate printed. |

## Reference

- **Wolfram code runs from a file.** `wolframscript -file script.wls` reads one expression at a time. An inline `-code` string is parsed whole before `Needs` runs, so every package symbol lands in the `Global` context and stays unevaluated.
- **`--timeout` belongs to `wolframscript`.** It acts on that argument even after the script name: the whole run stops after that many seconds, with exit code 0. The walker's own flag is `--solve-timeout`.
- **`getExampleScenario[key]` alone is a factory.** A scenario needs `getExampleScenario[key, entries, exits]`; `--case` supplies the defaults from `$ExampleBoundaryDefaults` in `Scripts/SolutionCacheHelpers.wls`.
- **The walk calls the stages directly** (`makeSymbolicUnknowns`, `makeSystem`, `dnfReduceSystem`), so the session cache of `solveScenario` plays no part.
