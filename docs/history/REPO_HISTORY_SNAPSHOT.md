# Repository History Snapshot

**Generated:** 2026-04-30
**Repository:** `MFGraphs`

## Core facts

| Metric | Value |
|---|---:|
| First commit | `2020-05-06` |
| Latest commit in snapshot | `2026-04-30` |
| Total commits on `master` history | **864** |
| 2020 commits | 160 |
| 2021 commits | 187 |
| 2022 commits | 89 |
| 2023 commits | 1 |
| 2024 commits | 18 |
| 2026 commits | 409 |
| March 2026 commits | 105 |
| April 2026 commits | 304 |
| March 2026 merged pull requests | 28 |
| April 2026 merged pull requests | 93 |

## Boundary commits

| Position | Commit |
|---|---|
| First | `09cc195` — `First commit` |
| Latest | `c594e8d` — `Merge pull request #175 from ricardodelimaribeiro/chore/ship-20260430-182409` |

## Representative recent history (April 2026)

| Date | Commit | Subject |
|---|---|---|
| 2026-04-30 | `c594e8d` | Merge pull request #175 from `chore/ship-20260430-182409` |
| 2026-04-30 | `b5568da` | docs: solver design notes, benchmark history, repo restructure |
| 2026-04-30 | `0edaaa6` | Merge pull request #168 from `codex/solver-diagnostics-work` |
| 2026-04-30 | `4e9f598` | fix solver branch-state semantics |
| 2026-04-29 | `f746e6b` | Merge pull request #166 from `codex/rename-symbolic-unknowns` |
| 2026-04-29 | `fb1577d` | Rename symbolic unknown bundle API |
| 2026-04-29 | `fc572f6` | feat: implement orchestration layer, fix tests, add validation, and update gitignore |
| 2026-04-29 | `621a89c` | default `solveScenario` to `dnfReduceSystem` |
| 2026-04-28 | `500cfba` | Improve graphics plot helpers |
| 2026-04-27 | `23653ba` | Relax zero-flow edge equations (#160) |
| 2026-04-27 | `71be300` | docs(CLAUDE.md): add single-test command, all three solvers, and benchmarking section |
| 2026-04-26 | `4bc0aa2` | Merge pull request #147 from `docs/reducesystem-benchmark-history` |
| 2026-04-25 | `6532d51` | Merge pull request #144 from `fix/system-cleanup-and-critical-congestion` |
| 2026-04-24 | `ae31734` | refactor: typed scenario kernel and modular load order |
| 2026-04-24 | `54f02bc` | fix(System.wl): patch symbol leak, tighten sub-record patterns, remove `EdgeList` |
| 2026-04-23 | `182b427` | Refactor: modularize `mfgSystem` builders and archive legacy components |
| 2026-04-23 | `694290f` | Archive `DataToEquations` and harden core-only guards |
| 2026-04-23 | `c1e51eb` | Reduce MFGraphs loader to scenario core and archive examples data |
| 2026-04-23 | `52f5379` | Add `CycleScenario`, `GraphScenario`, `AMScenario` public constructors |
| 2026-04-23 | `7be9bba` | Add `GridScenario` public constructor for direct grid/chain scenario creation |
| 2026-04-23 | `de335b9` | feat(examples): add self-contained `ExampleScenarios.wl` with concrete numeric scenarios |
| 2026-04-23 | `5246bcf` | feat(scenario): extend Hamiltonian defaults and edge params |
| 2026-04-22 | `2b2a177` | Extract system kernel and route `DataToEquations` through it |
| 2026-04-22 | `010d2b0` | Refactor scenario topology lifecycle and consolidate symbolic unknowns |
| 2026-04-22 | `7f3ee72` | Add typed unknowns head and enforce numeric scenario boundaries |

## Current active module surface

As reflected in the current repository tree and workflow guidance, the active package surface now centers on:

- `scenarioTools.wl`
- `examples.wl`
- `unknownsTools.wl`
- `systemTools.wl`
- `solversTools.wl`
- `graphicsTools.wl`
- `orchestrationTools.wl`

The current active fast test suite includes:

- `scenario-kernel.mt`
- `symbolic-unknowns.mt`
- `reduce-system.mt`
- `scenario-consistency.mt`
- `graphicsTools.mt`
- `orchestration.mt`
