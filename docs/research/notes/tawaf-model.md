# Tawaf MFG model

A modelling note for the Tawaf scenario builder (`MFGraphs/Tawaf.wl`),
companion workbook (`MFGraphs/TawafWorkbook.wl`), and tests
(`MFGraphs/Tests/tawaf.mt`). This file documents *what the model represents*
and *why the construction is shaped the way it is*. The package files
themselves are the source of truth for *how* the construction works.

## What this models

The ṭawāf (طواف) ritual at the Masjid al-Ḥarām in Mecca: seven
counter-clockwise circumambulations of the Kaʿba, starting and ending at
the al-Ḥajar al-Aswad (Black Stone). Pilgrims aim to salute the Black
Stone on each pass. The MFG layer captures the macroscopic flow problem:
how does the dense crowd distribute itself across the matāf pavement when
each pilgrim is solving an individual cost-of-time problem and congestion
is shared with everyone else on the same physical stones at the same
moment.

The default canonical instance is `makeTawafScenario[7, 8, 1]` — 7 rounds,
8 angular stations per round, 1 radial layer.

## The central modelling idea: unroll, then couple

A single physical loop is awkward to encode as an MFG state. The value
function needs to know how far along the ritual each pilgrim is — a
pilgrim on her first circuit standing at the Black Stone is in a
different ritual state than one on her seventh circuit at the same
physical location, so they should not share value-function values.

Unrolling each round into its own copy of the ring solves the state
problem: every logical node carries `(round, position)` and the value
function varies independently across rounds.

But the *physical* pavement is shared. Flow from round 1 and round 5
occupies the same paving stones simultaneously, so congestion must
couple. `makeTawafSystem` performs that coupling by grouping every logical
flow `j[u, v]` that maps to the same `(physical-segment, direction)`
identifier and rewriting each logical flow as the sum of its
round-shifted siblings. Concretely, on a 2×3 ring the forward segment
"position 1 → 2" carries logical flows `j[1,2]` (round 1) and `j[4,5]`
(round 2); the rewrite replaces each with `j[1,2] + j[4,5]` wherever it
appears in the Hamilton-Jacobi block (`EqGeneral`) and the switching
optimality block (`AltOptCond`).

The rewrite is intentionally limited to those two blocks — see [scope of
coupling](#scope-of-coupling) below.

This unroll-then-couple pattern is the central modelling idea. Everything
else in `Tawaf.wl` is bookkeeping.

## Coordinate codec

Layer-major flattening (`Tawaf.wl:22-23`):

```
vertex(r, p, l) = (l - 1) * rounds * nodesPerRound + (r - 1) * nodesPerRound + p
```

Layers stack outermost, then rounds within each layer, then positions
within each round. `tawafEncode` / `tawafDecode` are the single source of
truth; `TawafWorkbook.wl:112` re-derives the same encoding for the helix
plot.

## Topology produced by `makeTawafScenario[rounds, nodesPerRound, layers]`

- **Tangential edges**: forward only (counter-clockwise). Within a layer,
  position `p` connects to `p + 1`; the last position of a round connects
  to position 1 of the next round. The last position of the *last* round
  has no outgoing tangential edge — that node is the exit.
- **Radial edges**: bidirectional between adjacent layers at the same
  `(round, position)`. Only emitted when `layers > 1`.
- **Entries**: one per layer, at `(round=1, position=1)`, flow
  `100 / layers` (uniform partition of total flow 100).
- **Exits**: one per layer, at `(round=rounds, position=nodesPerRound)`,
  exit cost 0.
- **Switching costs**: empty.

## The Black Stone potential

Edges *originating* at any position-1 node (any round, any layer) carry
`EdgeV = -5.0`; all other edges carry `EdgeV = 0.0`
(`Tawaf.wl:131-137`).

The sign is what carries the meaning. A negative `V` is attractive in the
Hamilton-Jacobi formulation: pilgrims pay less to traverse edges leaving
the Black Stone, biasing trajectories toward saluting it on every pass.
The magnitude `-5.0` is a tuning constant chosen to be visible in solved
trajectories without overwhelming the congestion term; it is not
calibrated against a measured pace or ritual rate.

If you change the magnitude, expect the solved value-function gradient
near position 1 to shift correspondingly. The `EdgeAlpha`/`EdgeG`
override mechanism could be used to add finer ritual structure (e.g.,
attraction at the Yemeni Corner / Rukn al-Yamānī as well), but the
default builder does not.

## Layers

`layers > 1` produces concentric copies of the ring connected by radial
edges in both directions. The natural physical interpretation is
concentric pavement bands around the Kaʿba — pilgrims close to the
Kaʿba walk a short tight ring; those farther out walk a longer loose
ring; the radial edges represent lane changes.

The builder does not commit to this interpretation. Layers are a
generalisation handle, not a calibrated structural feature; treat the
multi-layer cases (`makeTawafScenario[r, n, 2]`, etc.) as exploratory.

## Scope of coupling

`tawafCouplingRules` (`Tawaf.wl:73-83`) only rewrites two of the system
blocks:

| Block                         | Rewritten? | Reason                                                                                                   |
|-------------------------------|------------|----------------------------------------------------------------------------------------------------------|
| `EqGeneral`                   | yes        | Hamilton-Jacobi sees physical congestion.                                                                |
| `AltOptCond`                  | yes        | Switching optimality at each logical vertex sees physical congestion on the candidate outgoing segments. |
| `EqBalanceSplittingFlows`     | no         | Balance lives at each *logical* node — round 1 vs round 5 conserve mass independently.                   |
| `EqBalanceGatheringFlows`     | no         | Same.                                                                                                    |
| `IneqJs`, `IneqJts`           | no         | Each logical flow has its own non-negativity.                                                            |
| `IneqSwitchingByVertex`       | no         | Vertex-level switching margins are logical.                                                              |
| `AltFlows`                    | no         | Logical-direction complementarity (`j[a,b] · j[b,a] = 0`).                                              |
| `IneqExitValues`/`AltExitCond`| no         | Boundary blocks; auxiliary entry/exit flows are not coupled (verified by `Tests/tawaf.mt:103-120`).      |

This is the right scope for the standard ṭawāf semantics. If you need
physical-edge complementarity (e.g., an opposing-direction lane sharing
the same stones) or vertex-level switching margins that see physical
sums, extend `makeTawafSystem` rather than the validator —
`isValidSystemSolution` reads whatever the system stores.

## Solver and validation

`makeTawafSystem` returns an ordinary `mfgSystem` with two blocks
rewritten in place. The standard solver pipeline applies unchanged:

```mathematica
s   = makeTawafScenario[2, 3, 1];
sys = makeTawafSystem[s];
sol = solveScenario[s];                 (* dnfReduceSystem by default *)
isValidSystemSolution[sys, sol]         (* True *)
```

End-to-end coverage in `Tests/tawaf.mt:143-160`. The 2×3×1 case is the
smallest non-trivial coupled instance (one shared physical segment,
`j[1,2] ↔ j[4,5]`).

## Density extension (opt-in)

The base Tawaf builder couples **flows** on shared physical edges. Flow
alone does not characterise congestion in the standard stationary MFG
formulation — the Hamiltonian carries a density `m` as well. The opt-in
density extension (`makeTawafScenario[…, "Density" -> True]` →
`makeTawafSystem`) augments the system with a per-edge density family
and the corresponding consistency block.

### Math

Critical-congestion Hamiltonian with `g(m) = -1/(2m)` and a strictly
negative potential `V(x) < 0`. The condition `H = 0` yields the
algebraic relation between the density and the signed flow on each
edge:

```
m = -(j² + 1) / (2 V)        equivalently        j = ±√(-2 m V - 1)
```

where `j` is the signed net flow `j[a,b] - j[b,a]` along the directed
pair on an undirected edge `{a,b}`. The `±` is the direction of net
flow; AltFlows kills one of the two non-negative directional variables
at equilibrium.

### What `makeTawafSystem` adds when `Density -> True`

- A symbolic family `m[{a, b}]` — one per non-boundary undirected
  logical edge, keyed by `Sort[{a, b}]`. Lives under `Ms` in the
  Hamiltonian sub-record (`systemData[sys, "Ms"]`).
- An inequality block `IneqMs` enforcing `m > 0`.
- A consistency block `EqDensityFlow`: one equation
  `j_phys² + 2 V_{ab} m_phys + 1 == 0` per shared physical segment,
  where `j_phys` and `m_phys` are the cohort sums (the existing flow
  coupling and a parallel density coupling are both applied; siblings
  collapse into a single physical constraint after `DeleteDuplicates`).
- A baseline strictly-negative potential on every non-boundary edge
  (`"BaselinePotential"` option, default `-1.0`); the Black-Stone
  discount `-5.0` is added on top of the baseline. `V` must be `< 0`
  everywhere `m` is defined or the consistency relation has no real
  solution.

### Modelling fork (current implementation = additive)

Two defensible reads of "the cost is `j = ±√(-2 m V - 1)`":

1. **Additive (current).** `EqGeneral` keeps its linear flow cost
   `j[a,b] - j[b,a]`. `m` enters only through `EqDensityFlow`, so the
   density is a first-class quantity tied algebraically to the
   physical flow sum but does not directly substitute into the HJB.
2. **Substitutive (deferred).** `EqGeneral`'s `j[a,b] - j[b,a]` is
   replaced by the signed `√(-2 m_phys V - 1)` expression. Tighter
   coupling; significantly harder for `Reduce` because every HJB
   equation gains a radical.

The current code implements the additive form. Switching to
substitutive is a one-line rewrite atop this baseline.

### Solver implications

The consistency block introduces quadratic terms in `j` and a
multiplicative coupling between `j_phys²` and `m_phys`. The DNF
prefilter on the existing complementarity blocks still applies; per-
branch reduce gets harder but stays in the symbolic pipeline.
`isValidSystemSolution` validates whatever the system stores, so no
solver-side change is required.

### Test coverage

`Tests/tawaf.mt` asserts (at the 2×3×1 size): five `m` symbols
generated, five `m > 0` inequalities, three `EqDensityFlow` equations
(one per physical segment, after the cohort dedup), and the explicit
shared-segment equation
`(j[1,2] - j[2,1] + j[4,5] - j[5,4])² + 2 V (m[{1,2}] + m[{4,5}]) + 1 == 0`.
A separate test confirms that the density block is absent on the plain
(non-density) builder so the existing pipeline is not perturbed.

## Visualisation

`TawafWorkbook.wl` provides three views:

- **Physical / augmented network plots** (`rawNetworkPlot`,
  `richNetworkPlot`) — the standard MFGraphs surface, applied to the
  unrolled topology. Useful for inspecting the coupled flow rewrites.
- **Helix plot** (`tawafHelixPlot`, defined in the workbook) — 3D layout
  where angle = position, height = `(round - 1) · nodesPerRound +
  (position - 1)` (so equivalent positions across rounds stack
  vertically), and concentric helices = layers. Dashed gray edges link
  same-position nodes across rounds, making the unroll visually
  interpretable.

## Provenance

This is a synthetic test problem **inspired** by the ritual to exercise
shared-physical-edge congestion in the MFG-on-networks framework. The
code does not cite a specific Tawaf-modelling paper. None of the seven
references in `docs/research/papers/` (Achdou 2023, Al Saleh 2024,
Bakaryan 2025, Cacace 2017, Camilli & Marchi 2016, Camilli 2015, Gomes
2019) are about the ṭawāf specifically — the seven are general
first-order MFG-on-networks references that the Tawaf scenario provides
a structured exercise for.

If a published Tawaf model exists in the literature you have in mind,
add the citation here and adjust the parameter choices (`-5` Black Stone
potential, uniform `100/layers` entry flow, exit cost 0) to match.

## Pointers

- Public API: `makeTawafScenario`, `makeTawafSystem` —
  `MFGraphs/Tawaf.wl`.
- Tests: `MFGraphs/Tests/tawaf.mt` (11 assertions, including end-to-end
  solve + validate).
- Interactive workbook: `MFGraphs/TawafWorkbook.wl` (sections for 2×3×1,
  3×4×1, 2×3×2, and a structure-only canonical 7×8×1).
- API summary: `CLAUDE.md` § "Tawaf scenario builder".
