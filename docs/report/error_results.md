# Error Results

Validation of the 2D solvers against the analytical solutions available in
`lbm::analysis` (`CouetteSolution2D`, `PoiseuilleSolution2D`). Every run
extracts the streamwise velocity along the vertical centerline with
`functional::extract_dx_profile_along_y_center`, and the profile is normalized
by the reference velocity `u0.dx` stored in the file header, so the plotted
curves are comparable in magnitude and not only in shape.

## Poiseuille flow

Setup (`simulations/openmp/poiseuille_flow_2d_{bgk,trt}.cpp`): channel with
rigid bounce-back walls on the bottom and top faces, pressure-periodic
inlet/outlet on the left and right walls, `Re = 100`, 100 000 iterations. The
flow is not driven by a moving wall: the pressure drop is set from the target
centerline velocity `u_max = 0.1` as

```
pin = pout + (nx / ny^2) * 8 * nu * u_max
```

so the exact solution to compare against is the parabola peaking at `u_max`.

Coarse channel -- 129x129, `Re = 100`:

![Poiseuille profile, 129x129](../imgResults/poiseuille_comparison_129_129_100.png)

Elongated channel -- 400x200, `Re = 100`:

![Poiseuille profile, 400x200](../imgResults/poiseuille_comparison_400_200_100.png)

What the two plots show:

- Both operators recover the parabolic profile: the simulated curves are
  parabolas that vanish at the two walls, as bounce-back prescribes.
- **BGK and TRT are indistinguishable at this scale** -- the blue BGK curve
  lies entirely underneath the orange TRT one. On a straight channel at
  `Re = 100` the TRT free relaxation parameter buys nothing: the two operators
  differ in how they treat the odd (anti-symmetric) moments, and this flow has
  no boundary layer curvature to expose the difference. The gain expected from
  TRT is on the wall placement of curved/obstacle boundaries, not here.
- The peak of the simulated profile settles at about `0.95 u_max` against the
  `1.0 u_max` of the analytical solution, i.e. a systematic ~5% deficit at the
  centerline, while the two curves agree closely in the near-wall region
  (`y/H < 0.2` and `y/H > 0.8`).
- The deficit **does not shrink when the channel is refined** from 129 to 200
  nodes across the height: it is therefore not a discretization error of the
  bulk scheme but a systematic offset between the imposed pressure drop and the
  reference `u_max` used to normalize -- the effective channel height seen by
  bounce-back walls (which sit half a lattice spacing outside the last fluid
  node) differs from the `H = ny - 1` used both in `PoiseuilleSolution2D` and in
  the `pin` formula above.

Since `compute_error()` returns the **absolute** global error against the same
`PoiseuilleSolution2D`, the reported `L2` error inherits this offset: it is
dominated by the constant centerline gap rather than by the local accuracy of
the scheme.

## Qualitative check

The time evolution of the velocity norm for the same runs, and for the Couette
counterpart, is collected in [`flow_animations.md`](flow_animations.md).

## Reproducing the plots

Run the simulation (from `build/simulations/`, so that the hardcoded `out/...`
paths resolve), then feed the profile files -- the BGK one, the TRT one and the
sampled analytical solution -- to the profile tool. Each file carries its own
`%%profile <model> <n> <u_ref>` header, and the header name is what ends up in
the legend:

```bash
python scripts/visualize_profile.py \
    out/data_poiseuille_openmp_129_100_01_bgk.bin \
    out/data_poiseuille_openmp_129_100_01_trt.bin \
    out/data_poiseuille_exact_129.bin \
    -o docs/imgResults/poiseuille_comparison_129_129_100.png
```

> The 400x200 runs are not one of the configurations hardcoded in the mains:
> they were produced by editing the `Config<2>` list (`grid_size`,
> `out_frames`, `out_data`) before rebuilding.
