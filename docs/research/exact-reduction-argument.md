# Exact reduction argument for the critical structural model

This note states the mathematical model and the obligations required for an exact
solver to preserve its whole solution set. It separates elementary equivalence
arguments from trusted computer-algebra operations and from checks on particular
saved outputs. It is not a proof-assistant certificate, a novelty claim, or a
claim that every accepted Hamiltonian produces a continuous MFG solution.

The implementation reviewed is the exact-foundation work based on commit
`f00fb27d557c6d253268c94a7e699ed0ccb2224f`. The independent review snapshot and
diagnostics are in
[`Results/exact-review/20261004T071139Z-d8f0deca/`](../../Results/exact-review/20261004T071139Z-d8f0deca/).
This argument states required behavior; it does not, by itself, certify that all
implementation paths satisfy it. In particular, the snapshot's unresolved-`Solve`
branch deletion was reproduced in the review diagnostics. Later corrections must
be identified by their source snapshot and checked separately.

## 1. Model, orientation, and domain

Let the physical graph be finite and undirected. Each physical edge has canonical
length one. For an edge with orientation `a -> b`, define

\[
 f_{ab}=j[a,b],\qquad b_{ab}=j[b,a],\qquad J_{ab}=f_{ab}-b_{ab}.
\]

The two directional currents are real and nonnegative, and

\[
 f_{ab}=0\quad\text{or}\quad b_{ab}=0.
\]

`u[a,b]` denotes the value at endpoint **b** on the physical edge `{a,b}`;
`u[b,a]` is its value at endpoint **a**. With the coordinate increasing from `a`
to `b`, set `U(0)=u[b,a]` and `U(1)=u[a,b]`. The generated critical edge equation is

\[
 u[a,b]-u[b,a]+J_{ab}=0,
 \qquad\text{equivalently}\qquad J_{ab}=U(0)-U(1).
 \tag{1}
\]

It is imposed at every physical edge, including zero-current edges. It does not
depend on plotted geometric length. A model with physical length `L` would instead
need `U(0)-U(L)=L J` under the same critical transport convention; silently using
the current unit-edge equation would change that model.

For an admissible transition `r -> v -> w`, let

\[
 t_{rvw}=j[r,v,w]\ge0,\qquad
 s_{rvw}=u[w,v]+c_{rvw}(t_{rvw})-u[r,v].
\]

The constraints are `s_rvw>=0` and `t_rvw=0 or s_rvw=0`, together with reverse
transition complementarity where both reverse transitions exist. Splitting and
gathering equations require each directional edge current to equal the sum of
its incident transition currents in the appropriate direction. Consequently mass
is conserved at the original vertices.

An auxiliary entrance has its prescribed incoming current. Its reverse boundary
direction is absent. An auxiliary exit has nonnegative outgoing current `q_e` and
a boundary value variable `v_e` satisfying

\[
 v_e\le\phi_e,\qquad q_e=0\quad\text{or}\quad v_e=\phi_e.
 \tag{2}
\]

The code keeps one canonical value variable on each auxiliary boundary edge,
with the boundary cost folded into (2); it does not impose the physical unit-edge
law on that auxiliary edge. An unused exit may therefore retain a value interval.
Selecting one value from that interval is an instance, not the complete structural
solution family.

The model considered here has exact `Alpha==1` on every physical edge, exact real
coefficients, and exact nonnegative entrance data. For the elementary linear
elimination proof below, coefficients are rational constants. The current
constructor permits some broader numeric affine coefficients; extending the
implementation guarantee to them requires the same exact algebra and real-domain
obligations. Nonlinear switching costs require a different solver argument.
`Infinity` denotes a blocked transition: its flow is zero and its switching
inequality is vacuous for finite real values. It is not an ordinary coefficient
that can enter finite linear algebra.

All claims concern the **completed scenario**. Symmetrizing adjacency, adding
auxiliary topology, metric-closing supplied switching costs, or rationalizing
inexact values can change the original user's mathematical data. These operations
must be disclosed with source provenance. Exact arithmetic on their output does
not establish equivalence to a distinct pre-completion or inexact input model.

## 2. Connection to the stationary MFG equations

The relevant local primary source is
[Al Saleh, Bakaryan, Gomes, and Ribeiro (2024), *First-order mean-field games on networks and Wardrop equilibrium*](papers/Al%C2%A0Saleh%20et%20al.%20-%202024%20-%20First-order%20mean-field%20games%20on%20networks%20and%20Wardrop%20equilibrium.pdf),
DOI `10.4171/PM/2124`.

Example 3.1, printed page 10, equation (3.4), derives

\[
 \frac{J^2}{2\rho^{2-\alpha}}+V(x)=g(\rho),
 \qquad
 J\int_0^1\rho^{\alpha-1}\,dx=U(0)-U(1),
 \tag{3}
\]

under the positive-density assumption `rho>0`. At `alpha=1`, the second relation
is precisely (1), including at `J=0`. The source's equations (4.3)–(4.10), printed
pages 17–19, supply the directional and transition complementarity, switching
optimality, splitting/gathering, and entrance conditions. The exit condition on
printed page 19 is an inequality with equality under positive exit flow. The
auxiliary entrance/exit construction and terminal-cost convention are described
in Remarks 2.1–2.2 and Sections 4.1–4.2.

The package's default Hamiltonian is

\[
 H(p,\rho)=\frac{p^2}{2\rho}-1+\frac1\rho,
 \quad\alpha=1,\quad V=-1,\quad g(\rho)=-\frac1\rho.
 \tag{4}
\]

For every real structural edge current `J`, define

\[
 \rho_J=1+\frac{J^2}{2}>0,\qquad
 U_J(x)=U(0)-Jx,\qquad 0\le x\le1.
 \tag{5}
\]

These formulas give a direct exact reconstruction:

\[
 H(U_J',\rho_J)
 =\frac{J^2/2+1}{1+J^2/2}-1=0,
 \qquad
 -\rho_J D_pH(U_J',\rho_J)=-U_J'=J.
 \tag{6}
\]

The transport current is constant, so its spatial derivative vanishes. Equation
(1) verifies the other endpoint. At zero current, `rho=1` and `U` is constant;
no division by the current or singular zero-density convention is needed.
Conversely, any positive-density classical edge solution of this default
Hamiltonian and transport system has constant `J=-U'`, hence (1), and `H=0`
forces the density in (5). This is an edgewise correspondence, combined with the
explicit structural vertex/boundary conditions above.

This recovery gives a narrow model-fidelity argument for scenarios that actually
use (4). To apply it to a saved collection, verify each scenario's Hamiltonian
metadata, physical edge lengths, and boundary/switching model. It does not follow
solely from the condition `Alpha==1`.

For general accepted `V/G`, the first equation in (3) remains an additional
positive-density solvability condition. For example, `V=1` and `g=0` would require

\[
 \frac{J^2}{2\rho}+1=0,\qquad\rho>0,
\]

which is impossible for every real `J`, including zero. A structural current/value
solution can nevertheless exist because the structural builder ignores these
metadata. It must not be certified as a solution of that Hamiltonian PDE.

The reference's Theorem 6.2, printed page 28, assumes positive costs, strictly
consistent switching costs, and regularity of all vertices in the associated
road-traffic model. Dormant edges and zero switching costs in this collection
prevent using that theorem without checking its hypotheses. The direct recovery
(5)–(6) avoids relying on that theorem for edge solvability; it does not establish
a general Wardrop correspondence, stochastic/time-dependent model, normalization
of total population, or a theorem for arbitrary affine tolls.

## 3. Solution-set invariant and reconstruction

Let `x` be **all original real unknowns**, including currents and endpoint values.
Let `F_original(x)` be the conjunction of the twelve original generated blocks:
entrance, splitting, gathering, critical edge relations, both nonnegativity
families, switching and exit inequalities, and the four complementarity families.
Do not substitute a pruned formula for `F_original` when checking this invariant.

For a rule list `R={x_i -> r_i}`, write

\[
 E_R(x)=\bigwedge_i(x_i=r_i).
\]

The invariant after each exact elimination step is

\[
 F_{\rm original}(x)
 \quad\Longleftrightarrow\quad E_R(x)\land C(x),
 \qquad x\in\mathbb R^n.
 \tag{7}
\]

`R` is not discarded after substitution. It recovers the eliminated coordinates;
`C` carries every remaining equality, inequality, disjunction, and parameter
restriction. The rules must be acyclic or otherwise explicitly normalized into
a well-defined simultaneous real reconstruction. Overwriting a definition is
permitted only if its previous relation is still retained or proved redundant.

Equivalently, divide the unknowns into retained parameters `y` and eliminated
variables `z`. For normalized real reconstruction `z=h(y)` and residual domain
`D(y)`, (7) becomes

\[
 \{x\in\mathbb R^n:F_{\rm original}(x)\}
 =\{(h(y),y):y\in\mathbb R^k,\ D(y)\}.
 \tag{8}
\]

Both directions are required:

1. Every original solution has retained coordinates in `D` and is recovered by
   `h`; this prevents losing branches or free parameters.
2. Every `y` in `D` reconstructs a real original solution; this prevents adding
   spurious solutions or certifying a complex-valued reconstruction over a
   purportedly real domain.

A variable mentioned in a domain may be fixed by a residual equality, so counting
residual symbols is not a dimension or nonuniqueness proof. The empty domain
represents infeasibility; an unresolved domain calculation does not establish
emptiness or nonemptiness.

### 3.1 Initial rules and affine equality elimination

A gathering equation is already a definition of a directional edge current as
a sum of transition currents. Replacing that current throughout the remaining
formula while retaining the definition satisfies (7). Symmetric zero-cost
switching inequalities imply `v1>=v2` and `v2>=v1`, hence `v1=v2`; taking a
canonical value in each equality component is another valid seed elimination.

For exact rational affine equations `A x=b`, row operations on `[A | b]` preserve
their solution set. In row-reduced form, each pivot variable is a constant-affine
function of nonpivot variables. Keeping these definitions and substituting them
into **the entire residual formula** preserves (7) in both directions. The
remaining inequalities restrict the free variables; row reduction alone does
not remove them.

A contradictory row `0=c` with nonzero `c` must make the residual false. A zero
row is redundant. Nonlinear coefficient arrays must not be truncated to their
linear part. Symbolic parameter-dependent pivots need case distinctions and
nonzero-denominator conditions; ordinary generic row reduction is not a proof
for such inputs. Rational constant coefficients avoid this last complication.

New definitions are composed into previous right-hand sides. A propagation
iteration limit is safe if all unprocessed constraints remain in `C`; it can
limit simplification without changing the set (8).

## 4. Net-current and zero-sum propagation

### 4.1 Fixed net current, including the degenerate case

Suppose the original constraints provide

\[
 f\ge0,\quad b\ge0,\quad(f=0\lor b=0),
\]

and the already justified reconstruction implies `f-b=q` for a rational constant
`q`. Then

\[
 f=\max(q,0),\qquad b=\max(-q,0).
 \tag{9}
\]

If `q>0`, the possibility `f=0` would force `b=-q<0`, so complementarity forces
`b=0` and `f=q`. If `q<0`, the symmetric argument forces `f=0,b=-q`. If `q=0`,
the two flows are equal and complementarity forces both to zero. Conversely,
the values in (9) satisfy nonnegativity, complementarity, and difference `q` in
all three cases.

Thus adding (9) preserves (7). It does not assert anything about a transition
routing not implied by other constraints. A solver may conservatively skip
nonrational or unresolved `q` and leave its original equations untouched.
It may not infer a direction from geometry, a scenario name, a numerical sample,
or a small numerical tolerance.

### 4.2 A vanishing sum of nonnegative flows

Suppose an original conservation equality, after substitution of previously
proved zero variables, is

\[
 \sum_{i=1}^k a_i t_i=0,
 \qquad t_i\ge0,
 \tag{10}
\]

with zero constant term and every `a_i>0`. Every summand is nonnegative; if any
`t_i>0`, its summand is positive, contradicting the zero sum. Hence every `t_i=0`.
If every `a_i<0`, multiply by `-1` and apply the same proof. Conversely all-zero
variables satisfy (10). This proves equivalence under the stated nonnegativity
premises.

Mixed signs do not justify this inference: `t1-t2=0` allows any equal nonnegative
pair. A nonzero constant does not justify it: `t1+t2=1` permits nonzero flows.
A missing nonnegativity premise does not justify it: `t1+t2=0` permits `(1,-1)`.

The polynomial used in (10) must actually be an original enforced equality.
If a typed system stores both raw balance polynomials and `EqBalance*` blocks,
the implementation must establish their consistency or derive (10) from the
enforced blocks. Assuming that correspondence after independent user edits to
the records can invalidate the inference.

### 4.3 Composition with existing definitions

The inferred equations are conjoined with `C` and then eliminated under (7).
For example, if `R` contains `f -> t1+t2` and (9) proves `f=0`, the new residual
must include `t1+t2=0`. Simply replacing the old rule by `f -> 0` would erase
the conservation relation. The conjunction-first method preserves it and lets
(10) subsequently derive `t1=t2=0` when their signs are known.

Induction on propagation iterations proves preservation of (7): each added
equation is a consequence of the current equivalent representation, and each
subsequent elimination preserves equivalence. No fixed-point convergence claim
is necessary. Stopping early is safe if the unresolved original residuals remain.

## 5. Disjunctions and finalization

After complete branch expansion, let branch `i` represent

\[
 B_i(x)=E_{R_i}(x)\land C_i(x).
\]

The full returned set is their union, represented by `Or_i B_i`. If `K` contains
rules common to every branch, factoring them gives

\[
 \bigvee_i B_i
 \quad\Longleftrightarrow\quad
 E_K\land\bigvee_i\big(E_{R_i\setminus K}\land C_i\big).
 \tag{11}
\]

Shared symbolic rules are equalities and must also remain represented. Substitution
through common definitions is valid under `E_K`. It does not authorize dropping
a branch merely because it has no additional restrictions.

In particular, `True or Q` is `True`; deleting a `True` branch can narrow the
union. `False or Q` is `Q`; deleting a proved-false branch is safe. An empty union
is false, not true. Equality extraction must leave every unextracted relation in
the residual and recover all eliminated coordinates.

For exact affine equations with constant coefficients, a successful complete
`Solve` result has one affine parameter family. Reading that family does not
select one nonlinear root. This reasoning does not extend to nonlinear equations
or to an unevaluated backend result. If `Solve`, `Reduce`, or a feasibility query
times out, fails, or remains unevaluated, the corresponding branch is **unknown**.
The whole solver must return a distinguishable unresolved/timeout result, or
retain that original branch unchanged. It must not replace unknown by `False`,
silently discard it from a union, or describe the surviving subunion as complete.

Likewise, pruning requires a proof of infeasibility, and branch deduplication
requires equality of branch content or another proved equivalence. A hash can
index candidate duplicates; hash equality alone is not an exact logical proof.

## 6. What is proved, what is trusted, and what remains open

The implications (9) and (10), the elimination invariant (7), and the Boolean
factorization (11) are elementary mathematical facts under their stated
hypotheses. The default-Hamiltonian reconstruction (5)–(6) is a separate exact
model argument. They establish no new general mathematical theorem.

An implementation-level certificate additionally depends on Wolfram evaluating
exact expressions, row reduction, substitution, `Solve`, and `Reduce` correctly;
on the actual code retaining every premise and domain; and on rejecting or
preserving unresolved backend outcomes. Tests can exercise these obligations,
including degenerate net current, free transition flow, unrestricted branches,
and failures. They do not prove all source code correct.

For a particular returned family `S`, exact soundness asks whether

\[
 \exists x\in\mathbb R^n:\ S(x)\land\neg F_{\rm original}(x)
\]

is false, and checks that `S` is nonempty. Completeness separately asks whether

\[
 \exists x\in\mathbb R^n:\ F_{\rm original}(x)\land\neg S(x)
\]

is false. External family parameters must be existentially projected when forming
`S` for the completeness question; soundness must cover every permitted parameter
value and a real-valued reconstruction. An infeasible original system is a
separate certified empty-set result.

The exact-foundation report records independent completeness proofs for two saved
cases and timeouts for eight others. A timeout is an unproved result, not a
counterexample and not a completeness certificate. A complete pipeline proof that
all its implementation obligations hold could justify completeness without a
second end-to-end `Reduce`; the elementary reduction identities alone do not
provide that proof for the entire baseline/current implementation.

Similarly, a measured execution-time ratio is evidence of runtime on specified
inputs and source snapshots. Calling it acceleration of the same complete
mathematical task also requires equivalence of the final returned sets. Exact
soundness of two outputs alone allows them to be different proper subsets and
does not establish that premise. The reported timings, mathematical validity,
model fidelity, completeness status, and artifact-content integrity should remain
separate claims.
