# Reproducibility and release

This guide describes the software checks, computational evidence, and release process for the critical-congestion MFGraphs solver used by *An Exact Solver for Stationary Mean-Field Games on Networks*.

## Scope

The active MFGraphs solver supports the **critical-congestion regime** only: every edge must have `Alpha = 1`. Non-critical systems deliberately fail rather than being treated as supported computations. Archived modules are outside this release path.

## Requirements

- Mathematica or Wolfram Engine 12.0 or later, with `wolframscript` on `PATH`.
- A clean clone of this repository.
- Sufficient memory for symbolic reduction. The Jamarat solve is materially more demanding than the merge, fork, and grid examples; record the operating system, Wolfram version, processor, memory, and command invocation for any timing comparison.

## Clean-clone validation

From the repository root, run:

```bash
wolframscript -file Scripts/RunTests.wls fast
wolframscript -file Scripts/CheckCriticalSurfaceTests.wls
```

The first command runs the active regression suite. The second confirms that the active test surface does not depend on archived or non-critical solver symbols.

## Manuscript artifacts live with the manuscript

Scripts that emit the manuscript’s PDFs deliberately live in the companion paper project’s `scripts/` directory, beside the `figures/` and `executables/` folders they update. This keeps MFGraphs focused on reusable research software and keeps paper assets and their generators versioned together.

For a clean verification run, invoke the paper scripts with this exact MFGraphs checkout:

```bash
wolframscript -file scripts/RegenPaperFigures.wls --output-dir /tmp/paper-figures
wolframscript -file scripts/RegenPaperResults.wls --mfgraphs-root /path/to/MFGraphs --output-dir /tmp/paper-results
wolframscript -file scripts/RenderJamaratFigures.wls --mfgraphs-root /path/to/MFGraphs --output-dir /tmp/jamarat-figures
```

The manuscript scripts also accept `MFGRAPHS_ROOT=/path/to/MFGraphs` in place of `--mfgraphs-root`. Before running them for a submission, check out the exact release tag or commit that the manuscript cites.

## Reproducing computation and timing evidence

Use the package benchmark harness for the road-merge, road-fork, and Jamarat systems:

```bash
wolframscript -file Scripts/PaperBenchmark.wls merge
wolframscript -file Scripts/PaperBenchmark.wls fork
wolframscript -file Scripts/PaperBenchmark.wls jamarat
```

The first two modes include timing repetitions and a profile-helper smoke check. Jamarat is timed once because it is a long symbolic computation. The pruning comparison used for the paper is available through:

```bash
wolframscript -file Scripts/BaselinePruningBenchmark.wls
```

Treat timing results as machine-dependent. Validate structural quantities, solver validity, and generated artifacts against the submission materials; do not expect identical wall-clock times on different hardware or Wolfram versions.

## Submission-release gate

Before an article submission or revision, maintainers should:

1. Start from a clean clone and run the active test and critical-surface commands above.
2. Check out the selected MFGraphs commit and run the paper project’s generators against that same checkout.
3. Inspect the regenerated files and record the execution environment used for performance numbers.
4. Create an annotated Git tag at the exact commit used for the manuscript.
5. Publish that tag as a GitHub release and archive it with a permanent DOI.
6. Confirm the included MIT license is appropriate for the release and update `CITATION.cff` with the archival DOI before publication.

The current repository URL is useful for readers, but the versioned release and DOI are the archival reference for a submitted manuscript.
