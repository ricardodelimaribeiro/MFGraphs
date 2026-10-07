# Exact computational foundation

This report concerns the stationary **critical structural system actually generated
by MFGraphs**, not a general nonlinear or time-dependent mean field game solver.
The numerical evidence and command ledger are in
[`Results/exact-foundation/`](../../Results/exact-foundation/). The curated input
definitions are [`Scripts/CuratedExactScenarios.wls`](../../Scripts/CuratedExactScenarios.wls).

**Demonstrated result:** ten selected networks have exact solutions or exact
rule/domain families whose every member satisfies all original generated
constraints. Selected solve times range from 0.005 to 28.069 seconds. Independent
completeness is proved for the competing-exit chain and the Camilli small topology;
the other eight completeness checks timed out and remain unproved. The full
[results table](../../Results/exact-foundation/curated-results.md) includes separate
construction, solve, soundness, and completeness-check times.

## Model and scope of exactness

`makeScenario` validates and completes a network, symmetrizes its adjacency into
physical connections, completes auxiliary entrance/exit topology, and caches that
topology. Integer vertex labels and disjoint entry/exit vertex sets are required.
`makeSymbolicUnknowns` constructs directional edge currents `j[a,b]`, transition
currents `j[a,v,b]`, and endpoint values `u[a,b]`. `makeSystem` assembles the original
constraint blocks. `solveScenario` normally calls `dnfReduceSystem` and memoizes
the result; the research runner calls system solvers directly and never times that
memoized retrieval.

For a physical edge `{a,b}`, let `m = j[a,b]-j[b,a]`. The generated critical relation
is `u[a,b]-u[b,a]+m == 0`, including at zero flow. Directional and transition flows
are nonnegative. Opposite physical directions and applicable reverse transitions
are complementary. Splitting and gathering equations conserve each directional
current. At a transition `(r,v,w)`, slack is
`u[w,v] + c[r,v,w](j[r,v,w]) - u[r,v]`; slack is nonnegative, and either the transition
flow or slack vanishes. Entrance rates are prescribed. An auxiliary exit value is
at most its terminal cost, and equals it when exit flow is positive. The auxiliary
value of an unused exit can remain free in an interval.

Every active structural solver requires exact `Alpha == 1` globally and on each
edge. A non-critical scenario can be constructed, but its structural solver returns
a `Failure`; it is not a non-critical solution. `V`, `G`, `EdgeV`, and `EdgeG` are
retained metadata, not density equations consumed by these solvers. Consequently
the computations establish statements about currents and endpoint values. They
do not reconstruct or validate continuous edge densities or a general Hamiltonian
PDE solution. Physical edges have the unit critical traversal law; plotted geometric
lengths are not travel-cost coefficients.

The collection uses exact integer/rational nonnegative inflows, rational terminal
and switching costs, and no oracle or numeric fallback. Scalar real costs,
`Infinity` (blocked transitions), and affine flow-dependent costs are accepted by
system construction; nonlinear switching functions are now rejected before linear
elimination. The paper collection does not establish general performance or
uniqueness results for affine flow-dependent tolls. `makeScenario` applies a metric
closure to explicitly supplied finite numeric switching costs before filling absent
costs with zero. A scenario named “Inconsistent” therefore need not be infeasible:
the completed model, not the pre-closure input or name, is the mathematical input.

Exactness has two independent requirements: exact source values and a valid exact
solution representation. `exactBoundaryValue` historically rationalizes inexact
numeric boundary inputs; this creates an exact rational problem but does **not**
prove the original measured/rounded value exact. The new report examines original
entry/exit metadata as well as equations and refuses to certify those inputs as
exact. Numeric diagnostic residuals, `solutionBranchCostReport`'s `0.5`/`N`, the
experimental flow-first `LeastSquares[N[...]]`, and numeric oracle routines are not
used as exact proofs. Literature metadata containing decimal reported values is
not used as a model coefficient.

| Result category | Evidence required |
|---|---|
| Exact determined solution | All unknowns assigned exact real values; every original constraint holds exactly. This alone does not prove uniqueness. |
| Exact parametric family | Exact rules plus an explicitly saved nonempty real parameter domain; every original constraint holds for **all** members and all branches. A domain can be a singleton, so this representation label alone does not establish nonuniqueness. |
| Complete solution set | In addition, an independent exact reference proves there is no original solution outside the returned set. |
| Unresolved | A returned expression, unevaluated backend, unproved domain, or validation timeout. It is not counted as a solved network. |
| Numerical approximation | Any inexact source coefficient or returned value. Rationalization does not upgrade provenance. |
| Exact infeasible | Independent original-system reduction proves `False`. A solver's empty set claim alone is insufficient. |
| Timeout / failure | The relevant stage did not finish; neither category implies infeasibility. |

The older flow-first `solutionResultKind` can label a transition-flow family
`Rules`; it remains a diagnostic compatibility label. Paper tables use
`exactSolutionReport` instead. A `FindInstance` point proves at most existence;
Boolean solvers with `ReturnAll -> False` deliberately return one surviving
branch. Graph-distance heuristics, symmetry restrictions, and oracle-selected
active sets are not completeness references for this collection.

## Correctness changes and proof obligations

The original validator returned `True` for `{}` on a network requiring inflow 10,
because it accepted every block that was not concretely false. It now requires
every block to be proved, and checks a parametric result under its residual domain.
The independent `exactSolutionReport` additionally checks exact provenance,
nonemptiness, and completeness against original generated blocks. It uses all
twelve conservation, boundary, sign, edge-relation, and complementarity blocks;
it does not validate against a preprocessed/pruned substitute system.

Both DNF finalizers formerly deleted `True` disjuncts. For example, the family
`x == 1 || (x == 1 && y >= 0)` could be narrowed to `x == 1 && y >= 0`.
Unrestricted branches now dominate their union. Exact quantified regressions
exercise this branch-loss case. `findInstanceSystem` and the Boolean solver paths
now distinguish timeout from infeasibility; unresolved feasibility checks retain
their arms. Boolean component/disjunct timeouts formerly became `False` or were
silently removed from a purported union. They now return `$TimedOut`, even if some
other branches finished. The affine builder guard and coefficient-array guard prevent
nonlinear terms from being silently discarded by a linear elimination routine.

`linearNetReduceSystem` is an explicit opt-in; no scenario-name dispatch and no
change of default routing is involved. It adds only consequences of the original
system:

1. If `f,b >= 0`, `f == 0 || b == 0`, and the accumulated exact equalities already
   give `f-b = m` for a rational constant `m`, then `f=Max[m,0]` and `b=Max[-m,0]`.
   The implementation checks the original nonnegativity/complementarity premises.
2. In a conservation equality, a sum of nonnegative variables with strictly
   positive coefficients can vanish only if each variable vanishes (likewise
   with all coefficients strictly negative). This propagates a zero directional
   current to its incident transition flows.
3. The new equations are conjoined with existing constraints before elimination.
   Existing edge-to-transition sum rules are never overwritten: doing that would
   discard conservation equations. The iteration is finite and bounded; stopping
   early leaves remaining constraints for the original DNF engine.

These implications preserve the full solution set. They do not select a preferred
transition routing. Their elementary identities are proved by real quantifier
elimination in tests, and small generated networks are compared to an independent
original-system reference. This is a correctness argument for the transformation,
not a claim of a new mathematical theorem or a new general solver algorithm.

Soundness validation substitutes returned rules into each original block, reduces
the exact parameter domain over the reals, and proves that the domain conjoined
with the negation of each block is empty. Thus a multi-branch family is checked
universally, not through one sampled point. Completeness uses original equalities
with Wolfram `Reduce` (not the package's `RowReduce` implementation), then bounded
Boolean expansion and exact per-conjunction reduction when tractable. No reference
branch is skipped or counted as infeasible on timeout. A timed-out reference remains
an explicitly unproved completeness claim. This is algorithmic independence within
Wolfram, not verification by a second CAS or a proof assistant.
Externally named family parameters are existentially projected before the
completeness query; soundness still covers every allowed parameter value. The
validator explicitly rejects non-critical models as `UnsupportedModel`.

## Curated networks and exclusions

| Case | Scientific purpose and interpretation |
|---|---|
| `competing-exits` | Three-vertex chain with inflow 100 and exit costs 0 and 10; exposes an unused expensive exit and its value family. |
| `diamond` | Existing example 12: four vertices with five physical connections, branching/merging and an internal cross-connection. |
| `triangle` | Three-vertex cycle with two equal-cost exits; a dormant cross-edge makes zero-current behavior explicit. |
| `braess-split` | Eight-vertex parallel-route topology with canonical rational switching tolls. |
| `braess-congest` | Seven-vertex merge/split bottleneck and nonunique transition routing. These two named examples are not by themselves a controlled proof of Braess's paradox. |
| `camilli-simple` | Two entrances on the four-vertex, four-edge Figure 1 topology from Camilli, Carlini, and Marchi (2015). |
| `grid3` | Nine-vertex structured cyclic network; transition-flow nonuniqueness. |
| `grid4` | Sixteen-vertex structured showcase; substantial redundant branching despite determined net edge currents. |
| `jamarat` | Two entrances and three exits; the selected result uses two exits and leaves the third unused. Its apparent residual parameter is fixed by an equality. |
| `case23` | A harder two-entrance/two-exit case retained specifically to test the previous lexicographic branch-state rescue. |

The Camilli topology was checked visually against Figure 1 on printed page 4187
(PDF page 15), and Section 5.1 states four vertices and four edges. The repository's
edge list `{1,2},{2,3},{2,4},{3,4}` matches the drawing after relabeling. The source
models a time-dependent stochastic meeting problem and uses edge geometry and a
time-dependent arrival cost. This collection instead prescribes stationary inflows
50 and 50, terminal cost zero, and unit critical edges. It is a literature-derived
topology, **not a replication of the paper's numerical result**. See the
[local source PDF](papers/Camilli%20et%20al.%20-%202015%20-%20A%20model%20problem%20for%20Mean%20Field%20Games%20on%20networks.pdf).

Candidates excluded or deferred, with evidence retained rather than erased:

- `Camilli 2015 general`: the importer omits two visually ambiguous segments to
  match the paper's stated edge count. The small verified topology gives cleaner
  provenance; no claim is made that the larger import reproduces Figure 3 exactly.
- `Achdou 2023 junction`: a valid finite stationary analogy, but adds less structural
  variety here and does not implement the paper's finite-horizon relaxed equilibria.
- `New Braess`: its custom congestion metadata is not a supported structural cost
  law in the current kernel, so the named economic interpretation would be misleading.
- “Inconsistent” shortcut examples: switching-cost closure changes the supplied
  costs. Their names are not acceptable independent infeasibility evidence.
- `Grid0505` and larger, `HRF Scenario 1`, and case 21: saved artifacts/workbooks
  report timeouts (and some kernel failures). They were not launched as an inventory
  sweep. The [HardCases workbook](../../MFGraphs/HardCases.wl) records the historical
  failures; its hand-derived/warm-start case-21 point is not a fresh complete solve.
  The retained 4×4 grid and case 23 provide a focused harder tier.

## Artifact integrity and experiment design

The first new harness trial is preserved as an **incomplete diagnostic run**. Its
original `SameQ` check failed. The diagnosis reproduced a Wolfram 15.0.1 behavior:
an Association containing a Graph compares unequal after WXF serialization even
though `Normal[association]`, each extracted value, graph vertices/edges/options,
and the printed exact expression agree. The solution, original constraints,
parameter assumptions, and validation report agreed directly. A minimal scalar
Association did not show the discrepancy. This isolates a container/opaque-Graph
representation comparison issue; there was no observed mathematical-content loss
or symbol-context change. The retained evidence has these exact differing paths:

| Payload path | Direct `SameQ` | `Normal` comparison | Printed exact expression |
|---|---|---|---|
| `Scenario / 1 / Model` | False | True | Identical |
| `Scenario / 1 / Topology` | False | True | Identical |
| `System / 1` | False | True | Identical |
| Solution, validation, assumptions | True | Not needed | Identical |

The minimal reproduction gives `False` for the Graph-containing Association,
`True` for its `Normal` form, and `True` for the extracted Graph. A preliminary
Graph-replacement/whole-container comparison also remained `False`; that failed
diagnostic is retained. We do not claim to know Wolfram's private internal cause,
or that raw `SameQ` now passes. The evidence supports replacing that inappropriate
container-identity test with the explicit content contract below. The original
failed artifacts, raw comparison outcomes, input forms, and mathematical checks
are linked in the results ledger.

The replacement in `Scripts/ExactArtifactIntegrity.wls` compares every Association
key/value, graph vertices/edges/options, sparse-array dimensions/entries/default,
and every other expression using an inert expression tree. It retains symbol
contexts and all payload metadata. It neither uses a tolerance nor treats hash
equality as mathematical-content equality. Regression tests deliberately alter
an exact solution by `1/10^30`, remove a domain boundary, change a context, cost
metadata, or Hamiltonian function, and require rejection. Hashes subsequently
detect changed file bytes. `ArtifactIntegrity` and mathematical `Soundness` /
`Completeness` are separate output fields.

Each runner invocation creates a unique directory, holds an advisory lock against
another curated invocation, and runs one Wolfram kernel at a time. Every sample
has the same bounded two-vertex DNF warm-up, then separately timed construction,
fresh solving, and exact validation. The extra preprocessing diagnostic is timed
after the solve, so its work does not warm the measured solve; solve time includes
its own preprocessing. Repetitions rotate source/method order and use fresh kernels.
After a solve timeout, later repetitions of that same case/source/method are
explicitly logged as skipped. Validation timeouts do not suppress solve repetitions.
An outer process watchdog terminates that sample's process group and records a
distinct process-timeout/incomplete-artifact outcome.

Manifests record the commit, working-tree changes, source checksums, Wolfram version,
machine information, exact case specification, options, commands, and per-sample
status. Final runs also snapshot all used package and runner sources, including
untracked new files. A completion marker, reopened-content check, file checksums,
and expected-row count are required; process exit status alone is insufficient.

## Reproduction

From the repository root with the licensed Wolfram kernel available:

```bash
wolframscript -file Scripts/RunTests.wls fast
python3 Scripts/run_curated_exact.py --methods linear-net --timeout 10 --validation-timeout 10 --tag curated
python3 Scripts/run_curated_exact.py --cases jamarat --methods linear-net --timeout 60 --validation-timeout 10 --tag jamarat
```

The first curated command is bounded to ten definitions, not the registered
inventory; Jamarat times out at its 10-second solve limit. The second curated
command performs only its justified 60-second rerun. To give every selected case
the same larger cap in one invocation, use `--methods linear-net --timeout 60`.
This is a cap, not a predicted completion time. Use each saved manifest's exact
command and source snapshot to reproduce the historical measurements.

Inspect a saved artifact without solving:

```bash
wolframscript -file Scripts/InspectCuratedResult.wls Results/exact-foundation/20261003T105818Z-curated-final-9dea72ad/008-grid4-current-linear-net-r1/result.wxf
wolframscript -file Scripts/AuditCuratedResults.wls Results/exact-foundation/selection.json Results/exact-foundation/my-new-audit
```

The audit refuses to overwrite an existing output directory. It rechecks exact
soundness, exports readable mathematics, and leaves saved completeness outcomes
unchanged. The preserved successful audit is
[`20261003T1134-artifact-audit`](../../Results/exact-foundation/20261003T1134-artifact-audit/).

## Measured results and mathematical interpretation

The selected runs used Wolfram **15.0.1 for Mac OS X ARM (64-bit), July 2, 2026**,
macOS 15.8.1, Apple M3 Max, 16 reported processors, and 48 GiB RAM. The initial
commit was `f00fb27d557c6d253268c94a7e699ed0ccb2224f`; the measured changes are
uncommitted on `research/exact-foundation`. Manifests and source archives identify
the actual code, rather than attributing the changes to that commit. Pre-existing
block-ordering script/CSV changes were preserved, along with the original baseline.

All rows below use `linearNetReduceSystem`. Times are seconds from individual
selected measurements; the controlled comparisons below retain every repetition.
The physical graph size excludes auxiliary boundary vertices. Exact content
integrity and exact soundness are **proved/passed for every row**; neither implies
completeness.

| Case | Vertices / edges | Mathematical output | Solve | Soundness | Completeness check | Completeness |
|---|---:|---|---:|---:|---:|---|
| Competing exits | 3 / 2 | Exit-value interval | 0.029343 | 0.000583 | 0.037945 | Proved |
| Diamond | 4 / 5 | Determined point | 0.006167 | 0.002176 | 10.102574 | Timeout |
| Triangle, two exits | 3 / 3 | Determined point | 0.036039 | 0.000618 | 11.115708 | Timeout |
| Braess split | 8 / 8 | Determined point | 0.009202 | 0.001285 | 10.315763 | Timeout |
| Braess congest | 7 / 8 | Transition-flow interval | 0.212292 | 0.336578 | 10.026626 | Timeout |
| Camilli simple | 4 / 4 | Determined point | 0.005104 | 0.001941 | 1.079243 | Proved |
| 3×3 grid | 9 / 12 | Transition-flow interval | 0.022748 | 0.004748 | 10.005202 | Timeout |
| 4×4 grid | 16 / 24 | Transition-flow polytope | 0.084974 | 0.016007 | 10.005757 | Timeout |
| Jamarat | 9 / 11 | Point encoded with an equality residual | 28.069179 | 0.003785 | 10.001212 | Timeout |
| Case 23 | 6 / 8 | Transition-flow interval | 1.824702 | 0.304603 | 10.005145 | Timeout |

The nominal completeness cap is 10 seconds; interruption overhead accounts for
slightly larger observed times. Construction took 0.003–0.049 seconds, so it was
not the limiting stage. Total validation time and all three unknown-family counts
are in the [CSV](../../Results/exact-foundation/curated-results.csv).

The saved parameter domains give concrete interpretable families. For competing
exits, `0 <= u["auxExit3",3] <= 10`, while the expensive exit carries zero flow.
For the 3×3 grid, `0 <= j[4,5,8] <= 25`; fixed edge currents admit different
transition assignments. Braess congest reduces to
`1/2 <= j[3,4,5] <= 201/4` with `j[3,4,6] == 201/4 - j[3,4,5]` and five other
residual variables zero. Case 23 reduces to `10 <= j[2,3,5] <= 50` with the other
five residual symbols fixed. Jamarat's only residual is `u["auxExit9",9] == 0`,
so its returned set is a singleton; the validator's conservative
`ExactParametricFamily` label describes its rule/domain representation, not a
claim of freedom. `FreeVariables` counts symbols before domain equalities, not
the dimension of the solution set. Full domains, including the 4×4 grid's coupled
inequalities, are preserved in the [readable exact files](../../Results/exact-foundation/20261003T1134-artifact-audit/).

## Controlled performance findings

The default branch-state threshold in the inspected code is zero. A forced
lexicographic experiment sets it to infinity; `Block-Edge` selects the branch-state
path through its non-default ordering. These are different from ordinary default
DNF. The saved interleaved ordering CSV contained 457 rows over 37 cases; its older
prose findings included a retraction where an ordering option had been inert.
The present conclusions use the fresh focused comparisons, not that prose alone.

| Case / comparison | Preserved default DNF | Current default DNF | Forced lexicographic | Block-Edge | New exact reduction |
|---|---:|---:|---:|---:|---:|
| 3×3 grid, two-repetition median | 0.106367 | 0.089816 | Not rerun | Not rerun | 0.023475 |
| 4×4 grid, two-repetition median | 29.452212 | 29.386865 | Not rerun | Not rerun | 0.087293 |
| Case 23 | Timeout at 60 s | Timeout at 60 s | 47.772972 | Not rerun | 1.881734 |
| Braess congest, two-repetition median | Not rerun | 0.224643 | 1.149425 | 0.823867 | 0.222689 |
| Jamarat, one measurement | Not rerun | 28.000112 | Not rerun | Not rerun | 28.069179 |

The 4×4 grid improved by about **337×** relative to the preserved default median
on identical exact inputs, warm-up policy, and 60-second solve limits. Case 23
changed from a censored 60-second default solve to approximately 1.9 seconds;
no finite default timing or exact speedup is inferred from that timeout. Both
default case-23 retries were explicitly skipped after their first timeout.
Lexicographic branch-state solving again rescued this case, at 47.462 and
48.084 seconds. There is no measured improvement for Jamarat or Braess congest.

For Braess congest, **Block-Edge versus lexicographic branch-state** measures the
ordering effect within the forced branch-state approach: 0.824 versus 1.149
seconds, about 28% lower time for Block-Edge. **Block-Edge versus default DNF**
compares different solver paths: 0.824 versus 0.225 seconds, about 3.67× slower
for Block-Edge. The first result is not evidence to prefer Block-Edge over the
default. These two small repetitions support a focused observation, not a
statistical generalization to other topologies or machines.

Ordinary preprocessing left 21 variables / 38 disjunctions for the 3×3 grid,
52 / 80 for the 4×4 grid, and 23 / 34 for case 23. Linear preprocessing already
determined all physical net currents on the grids and five of twelve physical/
auxiliary net currents in case 23. Solving still enumerated redundant alternatives.
Net pinning alone did not remove the major bottleneck; propagating zero sums of
nonnegative transition flows did. The slow stage on these cases was symbolic
branch enumeration after construction, rather than graph building. With the new
solver, independent completeness checking is now the dominant unresolved stage
for most core examples. Jamarat retains a solve-stage bottleneck around 28 seconds.

The repository-required `BenchmarkSystemSolver.wls --case grid-2x3 --timeout 5`
tagged runs were also preserved in isolated snapshots. Before/after direct repeat
times were 8.366 / 6.954 ms; these single microbenchmarks are workflow checks, not
the basis for the main speedup claim. They never timed memoized `solveScenario`
retrieval. See [BENCHMARKS.md](../../BENCHMARKS.md) and the
[run ledger](../../Results/exact-foundation/README.md) for exact commands and all rows.

## Verification and responsible paper claims

The final **fast** suite completed with **387 passed, 0 failed, exit 0** across
15 files. The separate `full` suite (which includes archived tests and the extra
fictitious-play suite) was not run. A pre-existing `m::shdw` warning in `tawaf.mt`
also appears in the initial 352-pass baseline; it is not a new execution error.
Focused reruns passed exact-validation 20/20, exact-artifact 10/10,
Boolean-minimize 21/21, and DNF-reducer 13/13. The critical-surface guard passed.
API generation checked 84 usage signatures successfully. Full logs and exit-status
reconciliation are in [`checks/`](../../Results/exact-foundation/checks/).
The report-placeholder guard, local-link check, and `git diff --check` passed.
The three touched module TeX files compile successfully; their retained logs include
layout/first-pass bookmark warnings, and no publication-ready PDF is claimed.
Existing inventory smoke tests allow some bounded solve timeouts, so their passing
status is not promoted to an exactness certificate for networks outside this collection.

A paper can responsibly claim existence of the saved exact equilibria, validity
of every member of the displayed parameter families under their saved domains,
and the measured runtimes on the identified machine and code snapshot. It can
claim a complete solution set for the competing-exit chain and Camilli simple,
and hence uniqueness for the latter exact point within the generated model.
The explicit non-singleton intervals/polytope certify nonuniqueness of those
returned solutions even where the absence of additional families is unproved.

Remaining blockers are specific: independent completeness timed out for diamond,
triangle, Braess split, Braess congest, both grids, Jamarat, and case 23; no
cross-CAS or proof-assistant validation was performed; the generated critical
structural model does not establish a non-critical, continuous-density, or
time-dependent MFG solution; and the named Braess examples do not establish a
controlled paradox comparison. Negative-input infeasibility is independently
covered by a regression, but no scientifically selected positive-input network
here is claimed infeasible. Larger excluded candidates have not been rescued by
this work. No novelty claim or paper manuscript is supplied.
