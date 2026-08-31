# Flow Animations

Animations of the velocity magnitude, rendered by `scripts/py/visualize_frames.py`
from the raw frame stream: `jet` colormap, origin at the bottom left, colour
scale fixed to the maximum norm reached over the whole run, so brightness is
comparable frame to frame. One frame is emitted every `niters / nframes` steps.

## Lid-driven cavity, 3D

`D3Q19`, 200x200x200, `Re = 1000`, lid on the `z = nz-1` face
([`configs/lid_cavity_3d.toml`](../../configs/lid_cavity_3d.toml)):

![3D lid-driven cavity, Re = 1000](../imgResults/lid_cavity_3d.gif)

Also available as [`lid_cavity_3d.mp4`](../imgResults/lid_cavity_3d.mp4), which
is a tenth of the size of the GIF -- prefer it wherever the target renders
video.

There is no tabulated 3D counterpart to Ghia in `benchmarks/`, so this case is
checked qualitatively here (the primary vortex forming under the lid, the
corner recirculations) and quantitatively through the exported centerline.

## Pipe flow, 3D

`D3Q19` on the CUDA backend, a `CylindricalShell` rasterized inside a box with
pressure-periodic inlet and outlet
([`configs/pipe_config.toml`](../../configs/pipe_config.toml) and its
variants):

![3D pipe, Hagen-Poiseuille](../imgResults/pipe_poiseuille.gif)

The steady profile this run converges to, and the error it scores against
`HagenPoiseuilleSolution3D`, are in [`error_results.md`](error_results.md).

## Reproducing an animation

`-o` is what turns the interactive window into a saved file; the `.mp4` writer
needs `ffmpeg` on the `PATH`, while a `.gif` extension falls back to `pillow`:

```bash
python scripts/py/visualize_frames.py out/frames_lid_cavity_3d.bin --title "Sim: 3D lid-driven cavity" -o docs/imgResults/lid_cavity_3d.mp4
```

> The 2D Couette and Poiseuille animations that this page used to list were
> dropped along with their frame files. Re-render them from a fresh run if they
> are wanted back: the commands above are the same, only the input file and the
> title change.
