import numpy as np
import matplotlib.pyplot as plt
import sys
import subprocess
import os


config_file = "weak_scaling_config.txt"

config_lines = [
    "# Configuration file for weak scaling\n",
    "Nx 50\n",
    "Ny 50\n",
    "Re 100\n",
    "U_lid 0.1\n",
    "Nsteps 5000\n",
    "Nframes 100\n",
    "out_norm out/temp.txt\n",
    "out_bench out/temp_bench.txt\n",
    "Col BGK\n",
    " ",
    "Nx 100\n",
    "Ny 100\n",
    "Re 100\n",
    "U_lid 0.1\n",
    "Nsteps 7500\n",
    "Nframes 100\n",
    "out_norm out/temp.txt\n",
    "out_bench out/temp_bench.txt\n",
    "Col BGK\n",
    " ",
    "Nx 150\n",
    "Ny 150\n",
    "Re 500\n",
    "U_lid 0.1\n",
    "Nsteps 7500\n",
    "Nframes 100\n",
    "out_norm out/temp.txt\n",
    "out_bench out/temp_bench.txt\n",
    "Col BGK\n",
    " ",
    "Nx 200\n",
    "Ny 200\n",
    "Re 500\n",
    "U_lid 0.1\n",
    "Nsteps 10000\n",
    "Nframes 100\n",
    "out_norm out/temp.txt\n",
    "out_bench out/temp_bench.txt\n",
    "Col BGK\n",
    " ",
    "Nx 250\n",
    "Ny 250\n",
    "Re 1000\n",
    "U_lid 0.1\n",
    "Nsteps 20000\n",
    "Nframes 100\n",
    "out_norm out/temp.txt\n",
    "out_bench out/temp_bench.txt\n",
    "Col BGK\n",
    " ",
    "Nx 300\n",
    "Ny 300\n",
    "Re 1000\n",
    "U_lid 0.1\n",
    "Nsteps 25000\n",
    "Nframes 100\n",
    "out_norm out/temp.txt\n",
    "out_bench out/temp_bench.txt\n",
    "Col BGK\n",
    " ",
]

runtimes = []

def write_config_file():
    with open(config_file, "w") as f:
        f.writelines(config_lines)

def run(max_cores : int) -> ():
    my_env = os.environ.copy()
    # run strong scaling
    for i in range(2,max_cores + 1, max_cores/6):
        print(f"Running with {i} cores")
        my_env["OMP_NUM_THREADS"] = str(i)
        cmd = ["./build/lbm-2-lbm", config_file]
        result = subprocess.run(cmd, capture_output=True, text=True, env=my_env)
        runtimes.append( end - start )
        # Access the output
        print("Output:", result.stdout)
        print("Errors:", result.stderr)
        print("Exit Code:", result.returncode)
        print(f"time: {runtimes[i-1]}")


def read_file():
    

def create_graph(max_cores : int) -> ():
    plt.figure(figsize=(10, 6))

    # Fix: x_axis needs to match the length of runtimes
    x_axis = np.arange(1, max_cores + 1)
    
    # Plot measured data
    plt.plot(x_axis, runtimes, 'o-', label='Actual Execution Time')

    plt.xlabel("Number of Cores")
    plt.ylabel("Time [ms]")
    plt.title("Strong Scaling Test")
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend()
    plt.show()

 
if __name__ == '__main__':
    if(len(sys.argv) < 2):
        print(f"Usage {sys.argv[0]} max_num_threads")
        sys.exit()
    max_cores = int(sys.argv[1])
    write_config_file()
    run(max_cores)
    read_file()
    create_graph(max_cores)

