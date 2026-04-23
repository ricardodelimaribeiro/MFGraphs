# MFGraphs Archive

This directory contains modules intentionally removed from the active core-only load path.

## Archived modules

- `DataToEquations.wl` — solver-era network-to-equation compiler and related symbolic pipeline helpers.
- `DNFReduce.wl` — solver-era symbolic DNF reduction helper.
- `FictitiousPlayBackend.wl` — solver-era internal fictitious-play backend.
- `Graphics.wl` — solver-era graphics helpers (`NetworkGraphPlot`, `SolutionFlowPlot`, etc.).
- `SolveMFGDispatch.wl` — solver-era `SolveMFG` method-dispatch layer.
- `Solvers.wl` — solver-era critical congestion solver and validation stack.

## Notes

- Archived files are kept for future restoration work.
- `MFGraphs/MFGraphs.wl` does not load modules from this directory in the current core scenario phase.
