# Architecture {#architecture}

[TOC]

The library separates four concerns that a naive LBM code usually mixes:
**what the domain is** (`lbm::Lattice`), **how it is evolved**
(`lbm::SolverBase` and its backends), **how a link that leaves a fluid node is
resolved** (`lbm::Solid`), and **where the data goes**
(`lbm::DataObservable` / `lbm::IDataListener`).

## Layers

| Directory | Namespace | Role |
|---|---|---|
| `include/lbm-sim/core/` | `lbm`, `lbm::utils`, `lbm::types` | Grid indexing, vectors, points, discrete velocity sets. |
| `include/lbm-sim/boundaries/` | `lbm::Solid` | Domain and obstacle boundary conditions; link resolution. |
| `include/lbm-sim/collision-detection/` | `lbm::CollisionDetection` | Analytic shapes (circle, box, airfoil, cylindrical shell) and the solid mask built from them. |
| `include/lbm-sim/collision-operators/` | `lbm` | `CollisionParams` (relaxation) and `CollisionStrategy` (BGK / TRT kernels). |
| `include/lbm-sim/solver/` | `lbm` | `SolverBase` and the `OpenMPSolver` / `CUDASolver` backends. |
| `include/lbm-sim/data/` | `lbm` | Listeners: `AsyncBinaryWriter`, `VtkWriter`, and the observable mixin. |
| `include/lbm-sim/analysis/` | `lbm::analysis` | Analytical solutions, error norms, the Ghia benchmark reader. |
| `include/lbm-sim/config/` | `lbm::config` | TOML parsing into `SimulationConfig`. See @ref configuration. |
| `include/lbm/` | `lbm::formatting`, `lbm::logging` | Backend-agnostic helpers shared with the rest of the project (CSV writer, logging facade). |

## Compile-time parameters

Four things are fixed when a simulation type is instantiated, and none of them
costs a branch or a virtual call at run time:

```cpp
lbm::LBMSimulation<dim, VelocitySet, cm_t>
lbm::OpenMPSolver <dim, VelocitySet, cm_t>
```

- `dim` — 2 or 3. `LBMSimulation` `static_assert`s that it matches
  `VelocitySet::dim`.
- `VelocitySet` — `lbm::D2Q9`, `lbm::D3Q19` or `lbm::D3Q27`. A velocity set is
  a struct of `constexpr` arrays: `dir`, `wi`, `opp`, plus `dim` and `ndir`.
- `cm_t` — `lbm::CollisionModel`, which also selects the matching
  `lbm::CollisionParams<dim, cm_t>` specialisation, so a parameter set cannot
  be paired with the wrong operator.
- `backend_t` — `lbm::ExecutionBackend`, carried by the solver and *deduced*
  by `LBMSimulation::solve()` from the solver it is handed. One simulation
  type therefore accepts either an OpenMP or a CUDA solver.

## Ownership

```
LBMSimulation ── owns ──> Lattice          (grid, u, rho, mask, obstacles, BCs)
      │         owns ──> CollisionParams  (Re, nu, tau)
      │
      └── solve(solver) ──> Solver ── owns ──> f_from, f_to  (populations)
```

The distribution functions are deliberately **not** in `lbm::Lattice`. They
belong to the solver, which picks their layout (AoS/SoA, host or device) to
suit its backend; the lattice holds only what every backend and every
post-processing step has to agree on. That is also why `solve()` takes the
solver by reference and does not store it: the same lattice can be handed to a
second solver, bearing in mind it still carries the first run's fields.

Every member of `lbm::Lattice` except `u` and `rho` is `const`, so a lattice
can be constructed or moved but never reassigned in place.

## The time step

`lbm::OpenMPSolver::solve()` allocates two population buffers plus a
`float` scratch array for the velocity norms, initialises the equilibrium, and
then repeats `niters` times:

1. `update_stream_collide()` — a single fused pass that streams (resolving each
   link through `lbm::Solid`), computes the macroscopic moments, evaluates the
   equilibrium and collides;
2. `std::swap(ffrom, fto)` — populations are read from one buffer and written
   to the other, so there is no aliasing and no extra copy;
3. every `nskips = niters / nframes` steps, `write_norms()` pushes a frame of
   velocity norms to the listeners.

The macroscopic fields are materialised only on the steps that need them — a
frame step, or the last iteration — which is what the `store_macroscopic` flag
on `update_stream_collide()` controls.

`nskips` is computed once in the `lbm::SolverBase` constructor, which also
rejects `nframes > niters` with `std::invalid_argument`.

## The two backends

`lbm::CUDASolver` follows the same five steps, because the physics is
literally the same code: the boundary kernels in
`boundaries/boundary-conditions.hpp` and the collision kernels in
`collision-operators/collision-strategy.hpp` are marked `LBM_HD_FUNC` and
compile for host and device alike. What differs is only the machinery around
them.

| | OpenMP | CUDA |
|---|---|---|
| Iteration | one `omp parallel for` over the flat cell index | one thread per node, block shape from `BackendProperties<CUDA>` |
| Out-of-domain work | none: the loop bound is exact | the launch grid is rounded up, so a thread outside the domain returns immediately |
| Velocity-set tables | read from the `constexpr` arrays | read from `__constant__` memory, uploaded once by `lbm::cuda::upload_lattice_constants()` |
| Buffers | two `std::vector`, swapped | two `lbm::cuda::DeviceBuffer`, swapped |
| A frame step costs | a `memcpy` into the listener buffer | a stream synchronisation plus a device-to-host download |

The `LBM_HD_FUNC` rule is what makes this work, and it constrains what may be
written in shared code: no allocation, no exceptions, no standard containers,
no iostreams. `lbm::detail::direction()`, `weight()` and `opposite()` in
`backend/utils.hpp` exist precisely to hide the storage difference behind one
call.

@warning The CUDA path does not currently compile: `solver/cuda-solver.cuh`
includes `lbm-sim/cuda/structs.cuh` and `lbm-sim/cuda/utils.cuh`, a directory
that no longer exists — those files moved to `lbm-sim/backend/cuda/`. The
backend is off by default (`LBM_ENABLE_CUDA=OFF`), which is why the build
stays green.

## Boundary handling

The domain sides carry an `lbm::Solid::DomainBC<dim>`: a raw array of
`2 * dim` bytes — 4 in 2D, 6 in 3D — regardless of the resolution. Immersed
bodies are described separately by a per-node `lbm::types::solid_mask_t`
holding an obstacle id (`lbm::types::FLUID` for a fluid node) plus a side
table `lbm::Solid::ObstacleData<dim>` mapping id to
`{boundary condition, wall velocity}`.

For each link `(p, dir)` leaving a node, `resolve_link()` produces an
`lbm::Solid::LinkResolution<dim>`: the already-wrapped source coordinate, the
boundary condition to apply, and the owning obstacle id. The available
conditions are

| Constant | Meaning |
|---|---|
| `lbm::Solid::NONE` | plain streaming from `src` |
| `lbm::Solid::BB_RIGID_WALL` | halfway bounce-back off a stationary wall |
| `lbm::Solid::BB_MOVING_WALL` | bounce-back with an imposed wall velocity |
| `lbm::Solid::PERIODIC` | coordinate wrap |
| `lbm::Solid::PRESSURE_PERIODIC_INLET` | wrap, then rescale to the inlet pressure |
| `lbm::Solid::PRESSURE_PERIODIC_OUTLET` | wrap, then rescale to the outlet pressure |

Resolution order matters and is fixed: a link that wraps does so *first*, then
a non-periodic face may claim it, and a moving wall wins over a rigid one on a
shared corner. A single mask slot per node cannot express this, which is why
the domain faces and the obstacle mask are kept apart. See the commentary in
`boundaries/utils.hpp` for the full argument.

`lbm::Solid::assert_consistent_domain_bc()` is called from the `LBMSimulation`
constructor and catches the common setup mistake of an axis that wraps on one
face and walls on the other.

## Output path

`lbm::LBMSimulation` and `lbm::SolverBase` both derive from
`lbm::DataObservable`, and they are *two distinct observables*: the simulation
emits the grid header once, the solver emits the frames. A writer that should
see the complete stream must therefore be attached to both.

```cpp
auto writer = std::make_shared<lbm::AsyncBinaryWriter>(path);
simulation.attachListener(writer);
solver.attachListener(writer);
```

`lbm::DataObservable` keeps non-owning handles: whoever creates the listener
stays responsible for its lifetime. Attaching and detaching are not
thread-safe with respect to dispatch — do both outside the run.

The byte layout every listener receives is documented in @ref output_formats.
