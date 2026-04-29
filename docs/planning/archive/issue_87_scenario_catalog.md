## Summary
After the scenario kernel is defined, MFGraphs should add a **scenario catalog** that registers named scenarios and supports controlled migration away from the current example-definition style.

The first migration wave should stay intentionally small and representative: one core benchmark case and one inconsistent-switching case first, with a paper-tier scenario only after the kernel assumptions are stable.

## Why this matters
- A catalog makes benchmark selection, documentation, and example discovery queryable by name, tier, tags, and eligibility.
- It provides a clean replacement for script-local case lists and scattered example metadata.
- It creates the bridge between the scenario kernel and a future benchmark v2 workflow.

## Scope
- define a `scenarioCatalog[]` or equivalent registry concept for named scenarios
- define `scenarioByName[...]` and the resolution rules for retrieving a canonical scenario object
- define a `deriveScenario[...]` planning contract for parent-child variants
- require lineage-aware metadata for derived scenarios, including parent identity and scenario-hash provenance
- migrate a minimal representative set first: a core 4x4-type case and an inconsistent-switching case; treat a paper-tier case as the next migration wave
- define a recursive merge policy for inheritance so nested blocks are not accidentally overwritten wholesale

## Deliverables
- a scenario catalog design note
- a migration checklist for representative scenarios
- a clear parent-child lineage policy for derived scenarios

## Out of scope
- bulk migration of the full examples library in one pass
- solver implementation changes beyond the catalog contract

## Related issues
- Depends on the scenario kernel issue.
- Enables a scenario-driven benchmark runner and documentation cleanup.
- Related to #82 and #83, but distinct from both.
