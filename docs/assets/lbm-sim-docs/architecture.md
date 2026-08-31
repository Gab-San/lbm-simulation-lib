# Architecture {#architecture}

This page describes how the pieces of `lbm-sim` fit together: who owns the
data, what each template parameter buys, what one time step actually does, and
where to cut in when adding a velocity set, a collision operator or a backend.

[TOC]

## The shape of the library

```
              simulations/*.cpp        one main = geometry + parameters
                     |
   CollisionDetection|                 LBMSimulation<dim, VelocitySet, cm>
   Shape -> CollisionArea  ------->      owns Lattice<dim> + CollisionParams
        |                                  |
        | Solid::compute_boundary_mask     | solve(solver, problem)
        v                                  v
   boundary_mask                    SolverBase<dim, VelocitySet, cm, backend>
   (1 byte per node)                  |-- MPISolver2D    (OpenMP)
                                      `-- CUDASolver2D   (CUDA)
                                             |
                                             | CollisionStrategy<dim, VS, cm>
                                             | Solid::apply_boundary_condition
                                             v
   DataObservable --notifyListeners--> IDataListener --> AsyncBinaryWriter
        ^                                                        |
        |-- LBMSimulation  (grid header)                         v
        `-- SolverBase     (velocity-norm frames)         out/*.bin  -->
                                                          scripts/visualize_*.py
```

Three responsibilities are kept strictly apart:

- **`LBMSimulation`** owns the state (`Lattice`) and the physics parameters
  (`CollisionParams`), and is the only class the user code talks to for
  post-processing (`compute_error()`, `compute_ghia_error()`, `output()`).
- **the solver** owns the *iteration*: the population buffers, the loop, the
  frame cadence. It borrows the lattice by reference for the duration of
  `solve()` and gives it back with the final macroscopic fields in place.
- **the listeners** own the *output*. Neither the simulation nor the solver
  knows what a file is; they emit `std::vector<char>` chunks.

## Compile-time parameters

Four parameters travel together through almost every template:

| Parameter | Values today | Enforced by |
|-----------|--------------|-------------|
| `dim` | `2`,`3` | `static_assert(dim == VelocitySet::dim)` in `LBMSimulation`; `static_assert`s in `Grid`, `CollisionParams`, `ErrorEvaluator`, `utils::ops::dot()` for `dim == 3` |
| `VelocitySet` | `D2Q9`,`D3Q19`,`D3Q27` | duck-typed: any struct exposing `dim`, `ndir`, `inv_cs2`, `wi[]`, `dir[]`, `opp[]` |
| `cm_t` (`CollisionModel`) | `BGK`, `TRT` | `if constexpr` chain in `CollisionStrategy::apply()`, with a `static_assert` on `MRT` |
| `backend_t` (`ExecutionBackend`) | `OPEN_MP`, `CUDA` | carried as a tag on `SolverBase`, so a solver of one backend cannot be passed where another is expected |

Nothing here is a run-time cost: `CollisionStrategy::apply()` resolves to a
single inlined kernel, and `LBMSimulation::solve()` is templated on the
backend tag so the solver call is a direct virtual dispatch on one known type
rather than a branch on the model.

The price is that an unsupported combination is a **compile error**, by design:
`LBMSimulation<3, D2Q9, BGK>` does not compile, and neither does a solver
instantiated with `CollisionModel::MRT`.

## Data model

### Grid and index maps

`lbm::Grid<dim>` is nothing but the extents plus the two maps used everywhere:

```cpp
scalar_index(p) = size.x * p.y + p.x;                    // one value per node
field_index(p, i, ndir) = ndir * (size.x * p.y + p.x) + i; // one value per link
```

`field_index` is **node-major**: the `ndir` populations of a node are
contiguous. That is what lets the inner loop pull a node into a local
`std::array<double, ndir>`, work on it in registers, and write it back once.

@note The same layout is used by the CUDA kernels, where it is a trade-off
rather than a win: consecutive threads handle consecutive nodes, so each read
of direction `i` is strided by `ndir * sizeof(double)` (72 bytes for `D2Q9`)
and global loads are not coalesced. The usual remedy is a direction-major
(structure-of-arrays) layout on the device; changing it means changing only
`Grid::field_index` and the kernels that call it.

`Grid::contains()` is the domain test used to decide whether a link is
interior or crosses the boundary.

### Lattice

`lbm::Lattice<dim>` holds the state that outlives a single time step:

| Member | Meaning |
|--------|---------|
| `grid` | extents and index maps |
| `u` | macroscopic velocity, one `Vector<double, dim>` per node |
| `rho` | macroscopic density, initialized to `rho0 = 1.0` |
| `boundary_mask` | one `boundary_t` (1 byte) per node -- `Solid::NONE` for fluid |
| `pin`, `pout` | inlet/outlet pressures for the pressure-periodic conditions |

The populations `f_i` are **not** in the lattice: they belong to the solver,
which allocates two buffers of `area * ndir` doubles and swaps them each step.
The lattice therefore stays the size of the macroscopic state, and the
expensive arrays live and die with `solve()`.

`boundary_mask` is 1 byte per node -- about 1.4% of the population arrays for
`D2Q9` -- which is why per-node boundary description is affordable at any
resolution.

## One time step

`MPISolver2D::update_stream_collide()` is a single parallel pass over the
nodes. Per node, in order:

1. **Pre-streaming density.** `r_wall = sum_i ffrom[p, i]`, the density at the
   node *before* streaming. Only the moving-wall condition uses it, but it is
   needed before any population is overwritten.
2. **Streaming (pull) with boundary resolution.** For every direction `i`, the
   source node is `src = p - c_i`. If `grid.contains(src)`, the population is
   pulled from the neighbour; otherwise the link leaves the domain and
   `Solid::apply_boundary_condition()` fills `fp[i]` according to the mask
   value *at `p`*.
3. **Macroscopic moments.** `rho = sum_i f_i`, `u = (sum_i f_i c_i) / rho`.
4. **Store.** `rho` and `u` are written back to the lattice only on frame steps
   and on the last iteration -- except on pressure-periodic nodes, which need
   fresh values every step because the inlet/outlet condition reads the
   neighbour's `rho` and `u`.
5. **Collision.** `CollisionStrategy::apply()` relaxes the local `fp` array
   towards equilibrium in place.
6. **Write back** into the second buffer, `fto`.

After the loop the two buffers are swapped (`std::swap` on the vectors: a
pointer exchange), and on frame steps the velocity-norm field is pushed to the
listeners.

Reading from `ffrom` and writing to `fto` is what makes the node loop
embarrassingly parallel: no node writes where another reads, so the OpenMP
`parallel for` needs no synchronization at all, and the CUDA kernel needs no
`__syncthreads()`.

Populations are initialized once, before the loop, by `init_equilibrium()`:
`f_i = w_i rho (1 + 3 c_i.u + 4.5 (c_i.u)^2 - 1.5 u.u)` evaluated on the
lattice's initial `rho = 1`, `u = 0`.

## Boundary conditions

A boundary condition is a rule for a *link*, not for a node: it says what
enters `fp[i]` when the source of direction `i` lies outside the domain. All of
them live in `lbm::Solid` (`boundaries.hpp`) as free `LBM_HD_FUNC` functions,
so host and device share one implementation.

| Constant | Function | What it does |
|----------|----------|--------------|
| `BB_RIGID_WALL` | `apply_bb_rigid_wall()` | half-way bounce-back: `fp[i] = ffrom[p, opp(i)]` -- the population that was leaving comes back reversed |
| `BB_MOVING_WALL` | `apply_bb_moving_wall()` | bounce-back minus the wall momentum term `2 w_opp rho (c_opp . u0) / c_s^2`, with `u0 = params.init_vel` |
| `PERIODIC` | `apply_periodic()` | pulls from the wrapped source `(size + p - c_i) mod size` |
| `PRESSURE_PERIODIC_INLET` / `_OUTLET` | `apply_periodic_with_pressure_variation()` | wraps like the periodic case, then re-anchors the pressure: `fp[i] = feq_i(p_target, u_wrapped) + f_i(wrapped) - feq_i(rho_wrapped, u_wrapped)` |
| `NONE` | -- | fluid node; nothing is applied |

The pressure-periodic pair is what drives the Poiseuille channel: the inlet
carries `lattice.pin`, the outlet `lattice.pout`, and the difference imposes
the pressure gradient in place of a body force.

@warning `apply_boundary_condition()` is reached **only** when
`!grid.contains(src)`, i.e. for links that leave the domain. A mask value on a
node in the interior of the fluid is therefore never acted upon: obstacles that
do not lie on the domain border are rasterized but not yet enforced. The
`FIXME` in `compute_boundary_mask()` marks the same gap.

## From geometry to mask

`lbm::CollisionDetection` turns analytic shapes into the per-node mask:

```
Segment / Circle / Parallelogram         (shapes, in a std::variant)
        |  getPerimeter()                 Bresenham rasterization, cached
        v
CollisionArea                            a shape group + its origin, one
        |  getPerimeter()                 obstacle id per area
        v
Solid::compute_boundary_mask()           perimeter nodes  ->  boundary_t
```

Each `CollisionArea` in the vector passed to `compute_boundary_mask()` gets an
id equal to its position, and `obst_type_map[id]` says which boundary type its
perimeter nodes carry. The rasterization uses integer Bresenham
(`algorithms/brasenham_rasterisation`), and the perimeters are cached in the
shapes (`mutable std::vector`), so building the mask for a large domain costs
one pass over the perimeter, not over the area.

@note `Shape<dim, DerivedShape>` is written as a CRTP interface, but its
methods are never instantiated: `CollisionArea` stores shapes in a
`std::variant` and reaches them through `std::visit`, which calls the derived
methods directly. The base class documents the expected interface rather than
enforcing it.

## Collision operators

`CollisionParams<dim, cm_t>` derives, once per run, everything the kernel
needs from the three physical inputs (`Re`, grid size, reference velocity):

```
nu = init_vel.dx * num_cells.y / Re          (lattice viscosity)
```

BGK then uses a single relaxation time, TRT a pair:

| Model | Stored | Relation |
|-------|--------|----------|
| `BGK` | `tauinv = 2 / (6 nu + 1)`, `omtauinv = 1 - tauinv` | equivalent to `tau = 0.5 + 3 nu`, `tauinv = 1 / tau` |
| `TRT` | `tauPlus = 3 nu + 0.5`, `tauMinus = 0.5 + 0.25 / (tauPlus - 0.5)`, `s_plus = 1/tauPlus`, `s_minus = 1/tauMinus` | fixes the magic parameter `(tauPlus - 1/2)(tauMinus - 1/2) = 1/4` |

The BGK constructor also guards stability on the host: it throws when
`tau <= 0.5` and warns outside `[0.55, 1.2]`.

`CollisionStrategy::apply()` dispatches on `cm_t` with `if constexpr`:

- **BGK** relaxes every direction towards equilibrium independently:
  `f_i <- omtauinv * f_i + tauinv * feq_i`.
- **TRT** works on opposite pairs. The loop visits each pair once (it skips
  `i > opp(i)`), splits populations and equilibria into their symmetric and
  antisymmetric parts, and relaxes them with `s_plus` and `s_minus`
  respectively; the rest direction, where `i == opp(i)`, has no antisymmetric
  part and is relaxed with `s_plus` alone.

Both kernels read the velocity set through `lbm::detail::direction()`,
`weight()` and `opposite()` (`backend.hpp`), which resolve to the host
`constexpr` tables or to the CUDA `__constant__` mirrors depending on
`__CUDA_ARCH__`. That indirection is the reason one collision implementation
serves both backends.

## The CUDA backend

`CUDASolver2D` mirrors `MPISolver2D` step for step; what changes is only where
the memory lives:

- `cuda::upload_lattice_constants<2, D2Q9>()` copies `dir`, `wi` and `opp` into
  `__constant__` memory once per `solve()`.
- Populations, mask, macroscopic fields and the norm buffer are device
  allocations; the host lattice is uploaded once at the start.
- The node loop becomes a 2D grid of `16 x 16` blocks sized with `ceil_div`, so
  the domain is covered exactly and every thread re-checks `grid.contains(p)`.
- `CollisionStrategy` is passed **by value** as a kernel argument -- it is
  trivially copyable by design, which is also why `CollisionParams` members are
  plain doubles and `LBM_HD_FUNC`-constructible.
- The buffer swap is a host-side pointer swap between kernel launches; all work
  is issued on one stream, and the stream is synchronized only on frame steps,
  when `rho`, `u` and the norms are copied back and handed to the listeners.

Everything else -- boundary functions, collision kernels, index maps -- is the
same code compiled for the device, thanks to `LBM_HD_FUNC`.

## Output pipeline and file formats

`DataObservable` is a mixin, not a base class with behaviour: both
`LBMSimulation` and `SolverBase` inherit it, so both can emit. Listeners are
non-owning `shared_ptr`s, and `notifyListeners()` copies the buffer for every
listener except the last, which receives it by move.

@warning `attachListener()` / `detachListener()` are **not** thread-safe with
respect to `notifyListeners()`. The contract is: attach before `solve()`,
detach after it returns. A writer attached to both the simulation and the
solver must be detached from both.

`AsyncBinaryWriter` is the only concrete listener: `acceptData()` pushes the
chunk onto a queue and returns, while a dedicated thread pops, writes and
flushes every 32 chunks. The destructor sets the stop flag and joins, so the
queue is fully drained before the file closes -- disk I/O overlaps with the
solve instead of stalling it.

Two file formats come out of a run, and both are consumed by the Python tools
in `scripts/`:

**Frames** (`AsyncBinaryWriter`, read by `scripts/visualize_frames.py`):

```
int32 nx | int32 ny                       written once by LBMSimulation
float32 norm[nx*ny]                       one block per frame, appended
```

**Centerline profile** (`LBMSimulation::output()`, read by
`scripts/visualize_profile.py`):

```
"%%profile <model> <n> <u_ref>\n"         ASCII header line
float64 profile[n]                        raw binary payload
```

`<model>` is `collision_model_to_string(cm_t)` and becomes the curve label in
the plot; `<u_ref>` is `params.init_vel.dx`, and the plotting script divides by
it, so profiles from different runs are directly comparable in magnitude. The
profile itself is produced by one of the extractors in `lbm::functional`,
passed in as a `std::function`.

## Analysis

Two independent facilities, both reached from `LBMSimulation`.

**Against an analytical solution** -- modelled on `dealii::VectorTools`:

```cpp
const analysis::PoiseuilleSolution2D exact(H, u_max);
const double err = simulation.compute_error(analysis::NormType::L2, exact);
```

`ErrorEvaluator<dim>::integrate_difference()` produces one magnitude
`|u_sim - u_exact|` per node, and `compute_global_error()` reduces it in the
requested norm (`L1` sums, `L2` sums the squares and takes the root,
`L2_squared` skips the root, `Linfty` takes the maximum). The result is the
**absolute** error: no division by the norm of the reference. Any new exact
solution is a subclass of `analysis::Function<dim>` overriding `value()`.

**Against the Ghia benchmark** (`analysis/benchmarks.hpp`, 2D only):
`compute_ghia_error()` reads a table, extracts the matching centerline from
`lattice.u`, normalizes it by the lid velocity, interpolates it linearly at the
17 tabulated coordinates and reduces the differences with the same
`ErrorEvaluator`. It returns a `NormErrorResult` carrying both the absolute
error and the relative one (divided by the norm of the reference values).

The table header must be `%%benchmark ghia <axis> <Re>`, where `<axis>` is `x`
or `y`; the parser compares `<Re>` against the simulation's Reynolds number and
throws on a mismatch.

@warning `benchmarks/ghia/bench_data_10000.txt` carries only
`%%benchmark ghia` -- two tokens -- so `parse_header()` indexes past the end of
its token vector. Use the `data_x_*` / `data_y_*` tables, or fix that header,
before calling `compute_ghia_error()` on it.

## Logging

`include/lbm/logging.hpp` wraps Quill: `setup_quill()` starts the backend
thread, `create_or_get_logger(name)` returns a named logger (`main`,
`simulation`, `solver`, `data_log`, `writer` are the names used across the
library). The header also registers `fmtquill` formatters and Quill codecs for
`Point` and `Vector`, which is what makes `LOG_INFO(logger, "{}", grid_size)`
work.

The active level is fixed at compile time by `LBM_LOG_LEVEL` (see
`cmake/Quill.cmake`): statements below it are compiled out, so a release run
pays nothing for the instrumentation left in the code.

## Extension points

| To add... | Do this |
|-----------|---------|
| **a velocity set** | write a struct exposing `dim`, `ndir`, `inv_cs2`, `wi[]`, `dir[]`, `opp[]` (`D3Q27` is already drafted, commented out, in `core/velocity-sets.hpp`). Note that `SolverBase2D` and both solvers currently hardcode `D2Q9`: generalizing them to the `VelocitySet` parameter is part of the work |
| **a collision operator** | add the enumerator to `CollisionModel`, specialize `CollisionParams<dim, cm>` with the relaxation rates, add the branch in `CollisionStrategy::apply()`, and extend `collision_model_to_string()` (it labels the profile files) |
| **a backend** | add an `ExecutionBackend` enumerator and derive from `SolverBase<dim, VS, cm, backend>`, implementing `solve(Lattice&, const CollisionParams&)`. Reuse `Solid::apply_boundary_condition()` and `CollisionStrategy` -- both are backend-agnostic |
| **an output format** | implement `IDataListener::acceptData()` and attach it to the simulation and the solver. Chunks arrive in emission order: the grid header first, then one velocity-norm block per frame |
| **an exact solution** | subclass `analysis::Function<dim>` and override `value()`; it plugs straight into `compute_error()` |
| **a shape** | add a class with `isCollidingWith()`, `contains()` and `getPerimeter()`, then list it in the `CollisionShapesT` variant |

## Known gaps

Collected here so they are not rediscovered from the source:

- **3D is scaffolding only.** `Grid<3>` index maps, `CollisionParams` for
  `dim == 3`, `ErrorEvaluator` and `utils::ops::dot()` all `static_assert`;
  `D3Q27` is commented out.
- **`problems/problem_3d.hpp` is stale.** It includes a header that does not
  exist in this tree (`lbm-2-lbm/problems/problem_base.hpp`) and refers to an
  `ExecutionBackend::MPI` that no longer exists. Nothing includes it, so it
  never reaches a compiler.
- **The problem abstraction is vestigial.** `LBMSimulation::solve()` takes a
  `const LidCavity2D&` by concrete type and never calls `init()` on it, so
  initial conditions other than "at rest" are not expressible yet.
- **Interior obstacles are not enforced** (see the warning above).
- **`MRT` is declared but not implemented.**
- **Solvers are fixed to `D2Q9`** through `SolverBase2D`, even though the
  boundary and collision code is already written against a generic
  `VelocitySet`.
