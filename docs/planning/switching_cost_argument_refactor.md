# Plan: switchingCosts argument refactor — list key to three-argument call

## Motivation

Currently switching costs are stored as an Association with 3-tuple list keys and looked up via `switchingCosts[{r, i, w}]`. The proposed change is to instead use three-argument calls `switchingCosts[r, i, w]`, making the API more consistent with how `j[r,i,w]` and `u[r,i]` are called, and eliminating the syntactic asymmetry between the switching-cost key and the triple it describes.

## Current representation

```wolfram
switchingCosts = <|{r1, i1, w1} -> cost1, {r2, i2, w2} -> cost2, ...|>
switchingCosts[{r, i, w}]   (* lookup *)
```

Keys are list triples; callers always wrap the triple in a list.

## Proposed representation

```wolfram
switchingCosts[r_, i_, w_] := cost   (* pattern-based definition *)
switchingCosts[r, i, w]              (* lookup — no list wrapping *)
```

This mirrors the natural reading: "switching cost from e_{r,i} to e_{i,w} at junction i".

## Affected sites

Every place that constructs, normalizes, or queries switching costs:

| File | Location | Change needed |
|---|---|---|
| `scenarioTools.wl` | `completeScenario`, `switchingCostsNumericQ`, normalization | Build definitions instead of Association entries |
| `systemTools.wl` | `isSwitchingCostConsistent`, `ineqSwitch`, `altSwitch`, `makeSystem` | `switchingCosts[{r,i,w}]` → `switchingCosts[r,i,w]` |
| `primitives.wl` or new file | `switchingCosts` symbol declaration | Declare with usage string |
| `API_REFERENCE.md`, `CLAUDE.md` | All `switchingCosts[{…}]` examples | Update to three-argument form |
| Test files | Any test that constructs or inspects switching costs | Update key construction |

## Trade-offs

**For:**
- Uniform call syntax with `j[r,i,w]` and `u[r,i]`
- No accidental list-vs-triple confusion at call sites
- Natural WL idiom — pattern-defined functions rather than Association lookup

**Against:**
- Association lookup is O(log n); pattern dispatch on a symbol with many definitions is O(n) in the number of defined triples (for large switching-cost tables this matters)
- `Normal[switchingCosts]` (to serialize or iterate all costs) no longer works; would need a separate registry or `DownValues[switchingCosts]`
- Requires clearing `switchingCosts` between scenarios (global state issue); the current Association-per-scenario approach is safer in multi-scenario sessions

## Recommended approach if implemented

Keep the internal storage as an Association (for O(log n) lookup and easy iteration), but expose a thin accessor `switchingCostAt[sc, r, i, w]` or overload `switchingCosts` as a curried function `switchingCosts[sc][r, i, w]` — so callers write natural three-argument syntax without the list wrapper, while the backing store remains an Association.

```wolfram
switchingCostAt[sc_Association][r_, i_, w_] := Lookup[sc, Key[{r, i, w}], 0]
```

This avoids global state, preserves O(log n) lookup, and removes the list-wrapping asymmetry.
