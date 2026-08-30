# Validation {#validation}

[TOC]

Four cases are used to check that the solver reproduces known physics. Three
have a closed-form solution and are scored against it directly; the fourth,
the lid-driven cavity, has none and is scored against the Ghia et al. (1982)
tables.

| Case | Dimension | Reference | Entry point |
|---|---|---|---|
| Couette flow | 2D | `analysis::CouetteSolution2D` | `LBMSimulation::compute_error()` |
| Poiseuille channel | 2D | `analysis::PoiseuilleSolution2D` | `LBMSimulation::compute_error()` |
| Hagen-Poiseuille pipe | 3D | `analysis::HagenPoiseuilleSolution3D` | `LBMSimulation::compute_error()` |
| Lid-driven cavity | 2D | Ghia et al. (1982) tables | `LBMSimulation::compute_ghia_error()` |
| Lid-driven cavity | 3D | none in-tree | export the profile, compare externally |

## Against an analytical solution

`compute_error()` takes only the exact solution and the norm; the approximate
field and the grid come from the simulation itself.

```cpp
const double err = simulation.compute_error(
    analysis::NormType::L2,
    analysis::PoiseuilleSolution2D(/*H*/ ny - 1, /*Umax*/ u_ref));
```

Internally this is the two-step split of `analysis::ErrorEvaluator`:
`integrate_difference()` builds a per-cell error field, then
`compute_global_error()` reduces it in the requested norm. Keeping the two
apart means the per-cell field can be exported to see *where* the error is,
and the norm can be changed without recomputing the difference.

The available norms are `L1`, `L2`, `L2_squared` and `Linfty`. The result is
**absolute and unnormalised**: there is no cell-volume weighting, because in
lattice units the cell measure is 1. Two different resolutions are therefore
only comparable after dividing by the cell count.

### Three things that make a converged run look wrong

**The wall is half a cell away from where you think.** With halfway
bounce-back the wall sits midway between the last fluid node and the first
solid node. For the pipe that means the effective radius is
`R = r_inner + 0.5`, where `r_inner` is what was handed to
`CylindricalShell`. Getting this wrong shifts the whole parabola.

**Solid nodes are included in the reduction.** `integrate_difference()`
iterates the entire grid, obstacle interiors included. Inside a wall the
solver leaves `u = 0`, so the exact solution has to return zero there as
well — which is exactly why `HagenPoiseuilleSolution3D` returns zero outside
the pipe instead of continuing the parabola into negative values.

**The transient is diffusive and long.** It scales as `L²/ν`. For the 129×129
Poiseuille channel that is `128²/0.129 ≈ 1.3e5` steps; below roughly `1e5`
iterations the profile has not settled and the L2 error stays high with
nothing actually wrong. The comments at the top of the shipped `.toml` files
carry the arithmetic for each case.

## Against the Ghia benchmark

```cpp
const analysis::NormErrorResult e =
    simulation.compute_ghia_error(LBM_BENCHMARKS_DIR "/ghia/data_y_100.txt",
                                  analysis::NormType::L2);
// e.relative, e.absolute, e.norm_type
```

The tables hold 17 points along a centreline per Reynolds number. The reader
checks the Reynolds number in the file header against `params.reyn_num` and
raises if they disagree, so a run cannot silently be scored against the wrong
table. The simulated centreline is interpolated linearly at the 17 tabulated
coordinates and divided by the lid velocity — Ghia's values are already
normalised that way, so it is the simulation that is rescaled, never the
reference.

Unlike `compute_error()`, this one returns **both** the absolute error and
the error relative to the norm of the reference, in a
`analysis::NormErrorResult`.

Defined for `dim == 2` only. Instantiating it on a 3D simulation is a compile
error; a 3D run that never calls it compiles fine.

## The 3D cavity

There is no analytical solution and no tabulated equivalent in this library.
A 3D cavity run is validated by exporting the centreline profile with
`functional::extract_dx_profile_along_z_center()` and comparing it against
external reference data:

```bash
python scripts/visualize_profile.py out/profile_cavity3d.dat --title "Re = 100"
```

## Checking the outputs structurally

```bash
python scripts/validate_outputs.py
```

Separate from the physics: it verifies headers, payload lengths, finite
values, and output paths declared twice by different sources. Run it before
trusting a plot.

## Reference

- U. Ghia, K. N. Ghia, C. T. Shin, *High-Re solutions for incompressible flow
  using the Navier-Stokes equations and a multigrid method*, Journal of
  Computational Physics 48 (1982) 387-411. A copy is in
  `docs/references/Ghia1982.pdf`, and the extracted tables are under
  `benchmarks/ghia/`.
