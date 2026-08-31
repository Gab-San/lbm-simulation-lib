# lbm-simulation-lib {#mainpage}

`lbm-simulation-lib` is a header-only C++17 framework for lattice Boltzmann
(LBM) fluid simulations. It runs the same simulation code on a multi-threaded
OpenMP backend and on an optional CUDA backend, in two or three dimensions.

At a glance:

| | |
|---|---|
| Dimensions | 2D and 3D |
| Velocity sets | `lbm::D2Q9`, `lbm::D3Q19`, `lbm::D3Q27` |
| Collision operators | `lbm::CollisionModel::BGK`, `lbm::CollisionModel::TRT` |
| Backends | `lbm::ExecutionBackend::OPEN_MP`, `lbm::ExecutionBackend::CUDA` |
| Boundary conditions | bounce-back (rigid and moving wall), periodic, pressure-driven periodic inlet/outlet |
| Distribution | header-only; `lbm::sim` is an `INTERFACE` CMake target |

The dimension, the velocity set, the collision operator and the backend are
all **compile-time** parameters, so no virtual dispatch happens inside the
time loop.

## Where to start

| Page | Content |
|---|---|
| @ref architecture | How the pieces fit together: lattice, solver, boundaries, collision, listeners. |
| @ref configuration | The TOML configuration format and how a main consumes it. |
| @ref output_formats | The frame stream, the profile file and the profiling CSV, byte by byte. |
| @ref validation | The four benchmark cases, the error norms, what makes a converged run look wrong, and the profiles the runs produced. |
| @ref performance | Timing a run, driving a scaling sweep, and the strong and weak scaling measured on a 64-core node. |
| @ref extending | Adding a velocity set, a collision operator, an obstacle shape or an output listener. |

The class and file reference generated from the headers is reachable from the
*Classes* and *Files* tabs. The narrative overview of the project and the
physical model live in the repository
[`README.md`](https://github.com/Gab-San/lbm-simulation-lib#readme); these
pages do not repeat them and focus on the API instead -- the measured results
are the exception, since the figures only mean something next to the code that
produced them.

## A minimal simulation

Every simulation main follows the same seven steps. The 2D lid-driven cavity,
stripped to its essentials:

```cpp
using namespace lbm;

// 1. Domain and boundary conditions: three rigid walls plus a moving lid.
const types::DimPoint<2> grid_size{129, 129};

Solid::DomainBC<2> dbc{};
dbc.low(0)  = Solid::BB_RIGID_WALL;   // x = 0
dbc.high(0) = Solid::BB_RIGID_WALL;   // x = nx-1
dbc.low(1)  = Solid::BB_RIGID_WALL;   // y = 0
dbc.high(1) = Solid::BB_MOVING_WALL;  // y = ny-1, the lid

// 2. No obstacle is immersed in the fluid: the mask is entirely types::FLUID.
types::solid_mask_t solid_mask = Solid::compute_solid_mask<2>({}, grid_size);

// 3. Collision parameters: nu and tau are derived from Re and the lid velocity.
const utils::Vector<double, 2> lid_velocity{0.1, 0.0};
CollisionParams<2, CollisionModel::BGK> params(100.0, grid_size, lid_velocity);

// 4. The simulation owns the lattice; the solver owns the populations.
LBMSimulation<2, D2Q9, CollisionModel::BGK> simulation(
    grid_size, std::move(solid_mask), {}, dbc, params);

OpenMPSolver<2, D2Q9, CollisionModel::BGK> solver(/*niters*/ 10000,
                                                  /*nframes*/ 100);

// 5. One writer, attached to both observables: the simulation emits the grid
//    header, the solver emits the frames.
auto writer = std::make_shared<AsyncBinaryWriter>("out/frames.bin");
simulation.attachListener(writer);
solver.attachListener(writer);

// 6. Run.
simulation.solve(solver);

// 7. Post-process: a centreline profile, and the error against Ghia et al.
simulation.output("out/profile.dat",
                  functional::extract_dx_profile_along_y_center);

const auto err = simulation.compute_ghia_error("benchmarks/ghia/data_y_100.txt",
                                               analysis::NormType::L2);
```

The complete versions of this and of the other cases are under
`simulations/openmp/` and `simulations/cuda/`.

## Building the documentation

Doxygen is picked up at configure time; the targets exist only when it is
found.

```bash
cmake -S . -B build
cmake --build build --target lbm-docs-full
```

- `lbm-sim-docs` — the library headers alone, in `docs/doxygen/lbm-sim/`;
- `lbm-docs-full` — this portal, in `docs/doxygen/`.

Both are configured by `docs/CMakeLists.txt`, which is added from the root
`CMakeLists.txt` when `LBM_BUILD_DOCS` is `ON` (the default for a top-level
build).

## Reference material

- `docs/references/Ghia1982.pdf` — the lid-driven cavity benchmark tables used
  by `lbm::analysis::compute_ghia_error()`.
- `docs/references/AMSC_LBM_hands_on_proposal.pdf` — the original assignment.
- `docs/report/` — the project report: measured errors, scaling, flow
  animations, the UML sketch (`LBM_UML.html`, editable as `LBM_UML.drawio`)
  and the theory notes.
- `docs/assets/` — interactive viewers for the
  [D3Q19](lbm_d3q19_directions.html) and
  [D3Q27](lbm_d3q27_directions.html) direction sets. The two links resolve in
  the generated site, where the viewers are copied next to the pages; reading
  this file in the repository, they are `docs/assets/*.html`.
