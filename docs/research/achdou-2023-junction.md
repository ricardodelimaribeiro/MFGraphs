# Achdou et al. 2023 Junction Note

Source: `docs/research/papers/Achdou et al. - 2023 - First order Mean Field Games on networks.pdf`

Achdou, Mannucci, Marchi, and Tchou study finite-horizon deterministic first-order mean field games on networks, focusing on a junction: several half-lines glued at one vertex. The MFGraphs registry entry `"Achdou 2023 junction"` is only a finite stationary analog of that topology. It uses one central vertex and four leaf vertices, with zero default switching costs and caller-supplied entries/exits.

## Repo-Relevant Takeaways

- Topology and movement policy should remain separate. The paper's junction defines incident edges at a shared vertex; movement choices at the vertex are policy/control data, not a reason to encode the graph as directed topology.
- Vertex behavior matters. Agents reaching the junction choose among incident edges, and future MFGraphs APIs may need an explicit representation for vertex-local choices rather than relying only on edge triples.
- Relaxed equilibria are trajectory-level objects. The paper's equilibrium is a probability measure over admissible trajectories, while the current scenario kernel works with stationary graph scenarios and symbolic edge/transition unknowns.
- The paper associates relaxed equilibria with mild solutions `(u,m)`, where `u` is a value function and `m(t)` is the time-dependent distribution induced by the trajectory measure.
- The continuity equation is treated weakly. This is relevant for future finite-horizon work, especially if mass can concentrate at vertices or propagate as singular measures.

## Current Scope Boundary

The imported example deliberately does not implement:

- Kakutani fixed-point machinery for relaxed equilibria.
- Probability spaces or measures over admissible trajectories.
- Viscosity-solution Hamilton-Jacobi theory on junctions.
- Finite-horizon first-order MFG solvers.
- Singular measure propagation or weak continuity-equation solvers.

## Future Design Candidates

These are possible names and concepts for later design work, not committed APIs:

- `VertexPolicy` for choices among incident edges and stay/leave behavior at vertices.
- Vertex and stay costs, separated from edge Hamiltonian data.
- Trajectory objects for finite-horizon dynamics.
- Relaxed-equilibrium records that carry trajectory measures and induced `(u,m)` mild solutions.
