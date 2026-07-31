import matplotlib.pyplot as plt
import sys
import subprocess
import os

# TO RUN WEAK SCALING BE SURE TO COMPILE 
# main.cpp WITH WEAK_SCALING to true

config_file = "weak_scaling_config.txt"

config_lines = [
    [
    "# Configuration file for weak scaling\n",
    "Nx 50\n",
    "Ny 50\n",
    "Re 100\n",
    "U_lid 0.1\n",
    "Nsteps 10000\n",
    "Nframes 100\n",
    "out_norm out/temp.txt\n",
    "out_bench out/temp_bench.txt\n",
    "Col BGK\n",
    ],
    [
    "Nx 100\n",
    "Ny 100\n",
    "Re 100\n",
    "U_lid 0.1\n",
    "Nsteps 10000\n",
    "Nframes 100\n",
    "out_norm out/temp.txt\n",
    "out_bench out/temp_bench.txt\n",
    "Col BGK\n",
    ],
    [
    "Nx 150\n",
    "Ny 150\n",
    "Re 500\n",
    "U_lid 0.1\n",
    "Nsteps 10000\n",
    "Nframes 100\n",
    "out_norm out/temp.txt\n",
    "out_bench out/temp_bench.txt\n",
    "Col BGK\n",
    ],
    [
    "Nx 200\n",
    "Ny 200\n",
    "Re 500\n",
    "U_lid 0.1\n",
    "Nsteps 10000\n",
    "Nframes 100\n",
    "out_norm out/temp.txt\n",
    "out_bench out/temp_bench.txt\n",
    "Col BGK\n\n",
    ],
    [
    "Nx 250\n",
    "Ny 250\n",
    "Re 1000\n",
    "U_lid 0.1\n",
    "Nsteps 10000\n",
    "Nframes 100\n",
    "out_norm out/temp.txt\n",
    "out_bench out/temp_bench.txt\n",
    "Col BGK\n\n",
    ],
    [
    "Nx 300\n",
    "Ny 300\n",
    "Re 1000\n",
    "U_lid 0.1\n",
    "Nsteps 10000\n",
    "Nframes 100\n",
    "out_norm out/temp.txt\n",
    "out_bench out/temp_bench.txt\n",
    "Col BGK\n",
    ]
]

runtimes = []
input_file = "weak_scaling_data.txt"

def run(max_cores : int) -> ():
    my_env = os.environ.copy()
    # run strong scaling
    for i in range(2, max_cores + 1, 2):
        with open(config_file, "w") as f:
            f.writelines(config_lines[i//2 - 1])
        print(f"Running with {i} cores")
        my_env["OMP_NUM_THREADS"] = str(i)
        cmd = ["./build/lbm-2-lbm", config_file]
        result = subprocess.run(cmd, capture_output=True, text=True, env=my_env)
        # Access the output
        print("Output:", result.stdout)
        print("Errors:", result.stderr)
        print("Exit Code:", result.returncode)


def read_file():
    with open(input_file, "r") as f:
        i = 0
        for line in f:
            runtimes.append(float(line.strip()));
            i += 1

def create_graph() -> ():
    # Fix: x_axis needs to match the length of runtimes
    x_axis = [2 * i for i in range(1, len(runtimes) + 1)]   
    # Plot measured data
    plt.plot(x_axis, runtimes, 'o-', label='Actual Execution Time')

    plt.xlabel("Number of Cores", fontsize=24)
    plt.ylabel("Time [s]", fontsize=24)
    plt.title("Weak Scaling", fontsize=24)
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend()
    plt.show()

 
if __name__ == '__main__':
    if(len(sys.argv) < 2):
        print(f"Usage {sys.argv[0]} max_num_threads")
        sys.exit()
    max_cores = int(sys.argv[1])
    run(max_cores)
    read_file()
    create_graph()
    # Destory the file
    if os.path.exists(config_file):
        os.remove(config_file)
        print("Config file removed.")
    else:
        print("Config file not found.")
    if os.path.exists(input_file):
        os.remove(input_file)
        print("Input file removed.")
    else:
        print("File not found.")
