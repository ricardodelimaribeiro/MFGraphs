# Scripts/repro/

Debugging and reproduction scripts preserved from active investigation sessions.
These are not part of the main package workflow — they are snapshots of one-off checks
that may be useful as starting points for future debugging.

## Subdirectories

| Subdirectory | Contents |
|---|---|
| `quick/` | Lightweight smoke checks and targeted feature verifications |
| `failures/` | Reproduction scripts for specific failure modes or regressions |

## Running scripts

All scripts assume they are run from the **repository root**:

```bash
# From the repo root:
wolframscript -file Scripts/repro/quick/quick_verify_inconsistent_examples.wls
wolframscript -file Scripts/repro/failures/repro_large_d2e_failure.wls

# Python helper (requires Results/inconsistent_examples_verification.json to exist):
python3 Scripts/repro/parse_inconsistent_verification.py
```

Scripts export outputs to `Results/` as they did when they lived at the root.
No script here is intended to be part of the test or benchmark suite.
