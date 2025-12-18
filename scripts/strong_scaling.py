import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import subprocess
import os


config_file = "strong_scaling_config.txt"

config_lines = [
    "# Configuration file for strong scaling\n",
    "Nx 129\n",
    "Ny 129\n",
    "Re 100\n",
    "U_lid 0.1\n",
    "Nsteps 10000\n",
    "Nframes 100\n",
    "out_norm out/temp.txt\n",
    "out_bench out/temp_bench.txt\n",
    "Col BGK\n"
]

runtimes = []
speedup = []

def write_config_file():
    with open(config_file, "w") as f:
        f.writelines(config_lines)

def run(max_cores : int) -> ():
    my_env = os.environ.copy()
    one_core_time = 0

    # TODO: calc speedup

    # run strong scaling
    for i in range(1,max_cores + 1):
        print(f"Running with {i} cores")
        my_env["OMP_NUM_THREADS"] = str(i)
        cmd = ["./build/lbm-2-lbm", config_file]
        start = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True, env=my_env)
        end = time.time()
        val = end - start
        runtimes.append( val )

        if(i == 1) {
            speedup.append(1.0);
            one_core_time = val 
        } else {
            speedup.append(one_core_time/val)
        }

        # Access the output
        print("Output:", result.stdout)
        print("Errors:", result.stderr)
        print("Exit Code:", result.returncode)
        print(f"time: {runtimes[i-1]}")


def create_graph(max_cores : int) -> ():
    # Fix: x_axis needs to match the length of runtimes
    x_axis = np.arange(1, max_cores + 1)

    fig, ax1 = plt.subplots()

    # --- First Axis: Execution Time ---
    color_time = 'tab:blue'
    ax1.set_xlabel('Number of Cores')
    ax1.set_ylabel('Execution Time (s)', color=color_time)
    line1 = ax1.plot(x_axis, runtimes, 'o-', color=color_time, label='Actual Execution Time')
    ax1.tick_params(axis='y', labelcolor=color_time)

    # --- Second Axis: Speedup ---
    # Create a twin axis that shares the same x-axis
    ax2 = ax1.twinx() 
    color_speed = 'tab:red'
    ax2.set_ylabel('Speedup', color=color_speed)
    line2 = ax2.plot(x_axis, speedup, 's-', color=color_speed, label='Speedup')
    ax2.tick_params(axis='y', labelcolor=color_speed)

    # --- Formatting ---
    plt.title("Strong Scaling Test: Time vs Speedup")
    ax1.grid(True, which="both", ls="-", alpha=0.2)

    # Combining legends from both axes
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc='upper center')

    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    if(len(sys.argv) < 2):
        print(f"Usage {sys.argv[0]} max_num_threads")
        sys.exit()
    max_cores = int(sys.argv[1])
    write_config_file()
    run(max_cores)
    create_graph(max_cores)
