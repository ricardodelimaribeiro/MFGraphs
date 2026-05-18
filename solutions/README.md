# Cached example solutions

One `.wxf` file per registered example scenario, holding the symbolic solution returned by `solveScenario` for the canonical boundary defaults declared in `MFGraphs/Tests/example-coverage.mt`.

Each file stores an `Association`:

```mathematica
<|
    "Key"       -> "Jamaratv9",
    "Entries"   -> {{1, 100}, {2, 100}},
    "Exits"     -> {{7, 0}, {8, 0}, {9, 0}},
    "Solution"  -> <| ... | ...rules... >,
    "Timestamp" -> "2026-05-12T11:48:29"
|>
```

## Refresh

```bash
wolframscript -file Scripts/RegenerateSolutions.wls               # all keys
wolframscript -file Scripts/RegenerateSolutions.wls Jamaratv9     # one key
wolframscript -file Scripts/RegenerateSolutions.wls --timeout 600 # custom per-scenario timeout
```

The regenerate script overwrites existing files. Solutions for scenarios that exceed the timeout are left absent (no file written).

## Test integration

`MFGraphs/Tests/example-coverage.mt` reads each cache file when present and asserts the fresh solve matches it; on a cache miss the test still asserts `scenarioQ` + `mfgSystemQ` and writes a fresh entry. Running the fast suite therefore populates missing caches as a side effect.

## Filename convention

- Numeric keys → `case_<n>.wxf` (`case_3.wxf`, `case_104.wxf`, ...).
- String keys → spaces and forward slashes replaced with `_` (`Jamaratv9.wxf`, `Camilli_2015_simple.wxf`, `Achdou_2023_junction.wxf`).
