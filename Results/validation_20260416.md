# Validation Note — 2026-04-16

This note records the **narrowly scoped** validation carried out on branch `manus/inconsistent-switching-examples` after the benchmark runner fix and the later decision to exclude `NonLinearSolver` and `Monotone` from the remaining checks.

## Scope used for this validation

The validation was intentionally limited to the smallest checks needed to confirm the current branch state.

| Area | What was checked | What was intentionally excluded |
|------|------------------|---------------------------------|
| Benchmark runner | Process-isolated execution with an explicit/resolved `wolframscript` path | Full-suite reruns |
| Example correctness | New inconsistent-switching examples and their expanded alternative-transition-flow clauses | Broad solver comparison |
| Paper-tier timing | `DataToEquations` and `CriticalCongestion` only, on `HRF Scenario 1` | `NonLinearSolver`, `Monotone` |

## Verified findings

### 1. Benchmark runner fix

The benchmark runner now resolves a robust process-kernel path and allows an explicit override to be threaded through the suite. This removed the earlier failure mode where process isolation reported a generic pipeline failure because the child `wolframscript` path was not resolved correctly.

A scoped large-tier process-isolated check completed successfully after this fix, confirming that the previous child-kernel launch failure was no longer blocking the large-case path.

### 2. Core case 9 regression investigation

The reported slowdown for **core case 9** was compared against an earlier reference commit using a temporary comparison worktree. The slower behavior was still present there, which indicates that the current branch's new inconsistent-switching example additions are **not** the cause of the observed timing change.

### 3. Inconsistent-switching examples

The new inconsistent-switching examples were checked against the `DataToEquations` alternative-transition-flow construction. The validation confirmed two things:

1. The examples are recognized as inconsistent with respect to switching-cost consistency.
2. Their generated alternative-transition-flow conditions contain extra cross-pair disjunctive structure beyond the standard consistent-cost reverse-pair pattern.

### 4. Narrow paper-tier check

A tightly scoped paper-tier run was executed only for `HRF Scenario 1`, and only for `DataToEquations` plus `CriticalCongestion`.

| Case | D2E time (s) | Critical time (s) | Critical result kind | Message |
|------|--------------|-------------------|----------------------|---------|
| HRF Scenario 1 | 0.07206 | 121.931556 | Failure | SymbolicTimeout |

The result artifact for this restricted run was written to `Results/paper_critical_only.json`.

## Workflow note

Two workflow constraints were applied deliberately during this validation:

1. **Baby steps only** — one minimal change at a time, followed by immediate verification.
2. **No full benchmark reruns unless truly necessary** — only narrow single-case or smoke checks were used.

## Current interpretation

The current branch is **unblocked at the benchmark-runner level**: process isolation can now launch correctly with an explicit or resolved kernel path. The remaining paper-tier bottleneck is now the actual symbolic workload of `CriticalCongestion` on `HRF Scenario 1`, not a runner failure.
