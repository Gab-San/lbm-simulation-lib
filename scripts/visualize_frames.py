import argparse
from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


def read_data(file_name):
    with open(file_name, "rb") as f:
        # reading grid size
        nx = np.fromfile(f, dtype=np.int32, count=1)[0]
        ny = np.fromfile(f, dtype=np.int32, count=1)[0]

        # reading velocities
        data = np.fromfile(f, dtype=np.float32)

    return nx, ny, data


def create_frames(nx, ny, data, num_iterations, vmax, output_file, save, title):
    fig, ax = plt.subplots()

    frame0 = data[0 : nx * ny].reshape(ny, nx)
    im = ax.imshow(frame0, cmap="jet", origin="lower", vmin=0, vmax=vmax)
    subtitle = ax.text(
        1, 1.015, "", ha="right", va="bottom", transform=ax.transAxes, fontsize=12
    )
    ax.set_title(title, loc="left")
    plt.colorbar(im, ax=ax, label="Velocity norms")

    def update(iter):
        frame_data = data[iter * nx * ny : (iter + 1) * nx * ny].reshape(ny, nx)
        im.set_data(frame_data)
        subtitle.set_text(f"Iteration n°{iter}")
        return [im, subtitle]

    anim = animation.FuncAnimation(
        fig, update, frames=num_iterations, interval=100, blit=False, repeat=True
    )

    plt.show()

    if save:
        writer = "ffmpeg" if str(output_file).endswith(".mp4") else "pillow"
        anim.save(output_file, writer=writer, fps=10)
        print(f"Plot saved in {output_file}")


def existing_file(path_str: str) -> Path:
    """Validate that the path exists and is a regular file."""
    path = Path(path_str)
    if not path.exists():
        raise argparse.ArgumentTypeError(f"File '{path_str}' does not exist.")
    if not path.is_file():
        raise argparse.ArgumentTypeError(
            f"Path '{path_str}' is a directory, not a file."
        )
    return path


class VerboseHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter
):
    pass


def parse_args():
    """Configure and parse command-line arguments"""

    description = """
    ===================================================================
                      FLUID SIMULATION VISUALIZATION TOOL
    ===================================================================
        """

    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=VerboseHelpFormatter,
    )

    parser.add_argument(
        "input_file",
        type=existing_file,
        metavar="FRAMES_DATA",
        help="Path to frames file.",
    )

    parser.add_argument(
        "--title",
        nargs="?",
        type=str,
        default="Sim: 2D Lid-Driven Cavity",
        help="Custom title",
    )

    parser.add_argument("-o", "--output", nargs="?", type=str, help="Output file path.")

    args = parser.parse_args()

    return args


def main():
    args = parse_args()

    nx, ny, data = read_data(args.input_file)

    # print("Stampo la norma delle velocità:\n", data)
    vmax = data.max()
    # il numero di iterazioni è guale alla lunghezza dei dati diviso nx*ny
    num_iterations = int(len(data) // (nx * ny))
    create_frames(
        nx,
        ny,
        data,
        num_iterations,
        vmax,
        args.output,
        save=args.output != None,
        title=args.title,
    )


if __name__ == "__main__":
    main()
