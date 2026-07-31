import sys

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


def create_frames(nx, ny, data, num_iterations, vmax, output_file, save=False):
    # create a figure and axis for the animation
    fig, ax = plt.subplots()
    ims = []

    for iter in range(num_iterations):
        # legge dati per ogni iterazione e salva in una matrice (ny, nx)
        frame_data = np.array(data[iter * nx * ny : (iter + 1) * nx * ny]).reshape(
            ny, nx
        )

        im = ax.imshow(frame_data, cmap="jet", origin="lower", vmin=0, vmax=vmax)
        title = ax.text(
            1,
            1.015,
            f"Iteration n°{iter}",
            ha="right",
            va="bottom",
            transform=ax.transAxes,
        )
        title.set_fontsize(12)

        ims.append([im, title])

    ax.set_title("Sim: 2D Lid-Driven Cavity", loc="left")

    # Add colorbar to show velocity norms with the label distanciated from the axis
    plt.colorbar(ims[0][0], ax=ax, label="Velocity norms")
    anim = animation.ArtistAnimation(fig, ims, interval=100, blit=False, repeat=True)
    plt.show()

    if save:
        # writer = animation.FFMpegWriter(fps=15, metadata=dict(artist='Me'), bitrate=1800)
        # anim.save('lid_driven_cavity_simulation.mp4', writer=writer)
        anim.save(output_file, writer="pillow")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage {sys.argv[0]} <input_norms> [<output_file_path>.gif]")
        sys.exit()

    input_file = sys.argv[1]
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    else:
        output_file = "lid_driven_cavity_simulation.gif"
    nx, ny, data = read_data(input_file)
    # print("Stampo la norma delle velocità:\n", data)
    vmax = max(data)
    # il numero di iterazioni è guale alla lunghezza dei dati diviso nx*ny
    num_iterations = len(data) // (nx * ny)
    create_frames(nx, ny, data, num_iterations, vmax, output_file, save=True)

