# Extending the library {#extending}

[TOC]

Four extension points come up in practice: a new **velocity set**, a new
**collision operator**, a new **obstacle shape**, and a new **output
listener**. None of them requires touching the solver.

## A new velocity set

A velocity set is a plain struct of `constexpr` members — there is no base
class and nothing to inherit. `lbm::D2Q9` is the shortest example. The
contract is:

| Member | Type | Meaning |
|---|---|---|
| `dim` | `lbm::types::dim_t` | spatial dimension; `static_assert`ed against `LBMSimulation`'s `dim` |
| `ndir` | `std::size_t` | number of discrete directions |
| `wi[ndir]` | `double` | lattice weights, must sum to 1 |
| `dir[ndir]` | `lbm::types::Direction<dim>` | the discrete velocities |
| `opp[ndir]` | `std::size_t` | index of the opposite direction of each `i` |

```cpp
struct D2Q5 {
  static inline constexpr lbm::types::dim_t dim = 2;
  static inline constexpr std::size_t ndir = 5;

  static inline constexpr double wi[ndir] = {2.0 / 6, 1.0 / 6, 1.0 / 6,
                                             1.0 / 6, 1.0 / 6};

  static inline constexpr lbm::types::Direction<dim> dir[ndir] = {
      {0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1}};

  static inline constexpr std::size_t opp[ndir] = {0, 3, 4, 1, 2};
};
```

`opp` is what bounce-back reads, so `dir[opp[i]] == -dir[i]` must hold for
every `i`; get it wrong and the walls leak rather than fail loudly. Index `0`
is the rest direction by convention throughout the library.

On the CUDA side the arrays are mirrored into `__constant__` memory by
`lbm::cuda::upload_lattice_constants<dim, VelocitySet>()`, which is
templated on the set and needs no change.

The direction numbering of the shipped 3D sets can be inspected interactively:
[D3Q19](lbm_d3q19_directions.html), [D3Q27](lbm_d3q27_directions.html).

## A new collision operator

Three pieces have to line up, all keyed on the same
`lbm::CollisionModel` enumerator:

1. **The enumerator** — add it to `lbm::CollisionModel` in `metadata.hpp`,
   *and* add its `case` to `lbm::collision_model_to_string()`. The trailing
   `std::to_string()` there is a fallback for out-of-range values, so a
   forgotten `case` stringifies as a number instead of failing to compile.
2. **The parameters** — specialise
   `lbm::CollisionParams<dim, cm_t>` with whatever relaxation rates the
   operator needs. The BGK specialisation stores `nu`, `tauinv` and
   `omtauinv`; the TRT one stores `tauPlus`, `tauMinus`, `s_plus`, `s_minus`.
   Because the specialisation is selected by `cm_t`, a parameter set can never
   be paired with the wrong operator.
3. **The kernel** — add an `apply_<model>()` to
   `lbm::CollisionStrategy` and dispatch to it from `apply()`. The
   `if constexpr` chain there ends in a `static_assert`, so an operator with
   parameters but no kernel is a compile error at the point of use rather than
   a silent no-op.

`apply()` is `LBM_HD_FUNC`, so the kernel is compiled for host and device
alike: keep it free of host-only facilities (no allocation, no exceptions, no
`std::` containers) if the CUDA backend is to support it.

The `MRT` enumerator exists but has no `CollisionParams` specialisation and no
kernel: it is a placeholder, and instantiating it fails to compile.

## A new obstacle shape

Shapes use CRTP, not virtual dispatch, so they can be stored by value in a
`std::variant` and visited without an indirect call. Derive from
`lbm::CollisionDetection::Shape<dim, Derived>` and provide two members:

```cpp
template <lbm::types::dim_t dim>
class Ellipse : public lbm::CollisionDetection::Shape<dim, Ellipse<dim>> {
  static_assert(dim == 2, "Ellipse currently supports only 2D");

  const lbm::types::Coordinate<dim> centre;
  const double a, b;   // semi-axes, in cells

public:
  Ellipse(lbm::types::Coordinate<dim> c, double a_, double b_)
      : centre(c), a(a_), b(b_) {}

  /// True when the lattice node is inside the solid.
  bool contains(const lbm::types::Coordinate<dim> &p) const {
    const double dx = (p.x - centre.x) / a;
    const double dy = (p.y - centre.y) / b;
    return dx * dx + dy * dy <= 1.0;
  }

  /// Inclusive integer bounding box, in the same frame as contains().
  lbm::CollisionDetection::AABB<dim> aabb() const {
    const int ha = static_cast<int>(std::ceil(a)) + 1;
    const int hb = static_cast<int>(std::ceil(b)) + 1;
    return {{centre.x - ha, centre.y - hb}, {centre.x + ha, centre.y + hb}};
  }
};
```

Then add the type to `lbm::CollisionDetection::CollisionShapesT<dim>`, the
`std::variant` that `CollisionArea` holds. Mind the dimension: `std::visit`
instantiates *every* alternative, so a 2D-only shape must not appear in the 3D
variant — that is why `Airfoil` and `CylindricalShell` sit on opposite
branches of the `std::conditional_t`.

Two conventions matter:

- `aabb()` is **inclusive on both ends and unclamped**;
  `lbm::Solid::compute_solid_mask()` clamps it against the grid. Returning a
  box that is too large only costs setup time — the mask is painted once — so
  err on the generous side. Returning one that is too small silently truncates
  the body.
- Shapes are defined **relative to the `CollisionArea` position**:
  `CollisionArea::contains()` calls `shape.contains(point - position)`.

An obstacle only becomes a wall once its id is given a boundary condition in
the `lbm::Solid::ObstacleData` table passed to `LBMSimulation`; ids are the
indices of the areas handed to `compute_solid_mask()`.

## A new output listener

Implement `lbm::IDataListener`, whose single method receives raw byte chunks:

```cpp
class MyListener : public lbm::IDataListener {
public:
  void acceptData(std::vector<char> data) override;
};
```

The chunks arrive in the order described in @ref output_formats: the grid
header first, then one frame per emission. A listener that needs the header
must therefore be attached to the simulation *and* to the solver, since these
are two distinct observables:

```cpp
auto listener = std::make_shared<MyListener>();
simulation.attachListener(listener);
solver.attachListener(listener);
// ... run ...
simulation.detachListener(listener);
solver.detachListener(listener);
```

Notes worth knowing before writing one:

- `lbm::DataObservable::notifyListeners()` copies the buffer for every
  listener except the last, which gets it moved. Listener order is attach
  order.
- Attach and detach are **not** thread-safe with respect to dispatch: do both
  outside the run.
- The observable keeps a non-owning handle, so the listener must outlive the
  run. `AsyncBinaryWriter` and `VtkWriter` both drain their queue in the
  destructor; a listener that buffers must do the same or lose the tail.
- `acceptData()` runs on the calling thread. If the work is not cheap, queue
  it — that is exactly what `AsyncBinaryWriter` does with its worker thread.

## Adding a simulation executable

`simulations/CMakeLists.txt` globs `*.cpp` (and `*.cu` when
`LBM_ENABLE_CUDA=ON`) recursively and creates one target per file, CUDA ones
prefixed `cuda_`. Dropping a new main under `simulations/openmp/` is enough;
it automatically gets `lbm::sim`, plus the `LBM_CONFIGS_DIR` and
`LBM_BENCHMARKS_DIR` defines.

The glob being recursive, two sources with the same file name in different
directories would collide. That is caught at configure time with an explicit
message naming both, rather than as a duplicate-target error further down.
