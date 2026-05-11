# MFGraphs Archive

This directory contains modules intentionally removed from the active core-only load path.

## Archived modules

- Solver-era modules and backend helpers that are not loaded by the active package.
- `getKirchhoffMatrix.wl` — public wrapper around `getKirchhoffLinearSystem` with a placeholder cost-function slot; archived after the active code path stopped consuming the 4-tuple shape.

## Notes

- Archived files are kept for future restoration work.
- `MFGraphs/MFGraphs.wl` does not load modules from this directory in the current core scenario phase.
