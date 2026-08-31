# Flow Animations

Animations of the velocity norm produced by the 2D OpenMP solvers, one frame
every `niters / nframes` steps (200 frames for the reference configuration:
100 000 iterations, `nframes = 200`). They are rendered by
`scripts/visualize_frames.py` from the raw `AsyncBinaryWriter` output: `jet`
colormap, origin at the bottom-left, colour scale fixed to the maximum velocity
norm reached over the whole run, so brightness is comparable frame to frame.

All runs are at `Re = 100` with reference velocity `u0 = (0.1, 0)`; `129x129`
is the square reference domain, `400x200` the elongated channel.

## Couette flow

Rigid bounce-back bottom wall, moving top wall, periodic left/right faces. The
animation shows the linear profile building up from rest as momentum diffuses
down from the lid.

| Grid | BGK | TRT |
|------|-----|-----|
| 129x129 | [couette_129_bgk.mp4](../imgResults/frames_couette_129_129_100_01_bgk.mp4) | [couette_129_trt.mp4](../imgResults/frames_couette_129_129_100_01_trt.mp4) |
| 400x200 | [couette_400x200_bgk.mp4](../imgResults/frames_couette_400_200_100_01_bgk.mp4) | [couette_400x200_trt.mp4](../imgResults/frames_couette_400_200_100_01_trt.mp4) |

## Poiseuille flow

`D2Q9`, rigid bounce-back top and bottom walls, pressure-periodic inlet and
outlet (see [`error_results.md`](error_results.md) for the pressure drop). The
animation shows the parabolic profile developing from the uniform initial
condition until the pressure gradient and the viscous stresses balance.

| Grid | BGK | TRT |
|------|-----|-----|
| 129x129 | [poiseuille_129_bgk.mp4](../imgResults/frames_poiseuille_d2q9_129_129_100_bgk.mp4) | [poiseuille_129_trt.mp4](../imgResults/frames_poiseuille_d2q9_129_129_100_trt.mp4) |
| 400x200 | [poiseuille_400x200_bgk.mp4](../imgResults/frames_poiseuille_d2q9_400_200_100_bgk.mp4) | [poiseuille_400x200_trt.mp4](../imgResults/frames_poiseuille_d2q9_400_200_100_trt.mp4) |

> The files are `.mp4`: GitHub plays them inline once the link is opened, and
> any local player handles them directly. As with the profiles, BGK and TRT are
> visually indistinguishable on these two flows.

## Reproducing an animation

`-o` is what turns the interactive window into a saved file; the `.mp4` writer
needs `ffmpeg` on the `PATH` (drop the extension to `.gif` to fall back to
`pillow`). Paths are the frame files written by the run, relative to
`build/simulations/`:

```bash
python scripts/visualize_frames.py out/norms_poiseuille_openmp_129_100_01_bgk.bin --title "Sim: 2D Poiseuille flow" -o docs/imgResults/frames_poiseuille_d2q9_129_129_100_bgk.mp4
```
