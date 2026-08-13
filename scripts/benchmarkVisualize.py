import sys

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


def read_data(file_name, bench_name):

    # leggo i dati di ghia
    with open(bench_name, "r") as f:
        ghia_data = []
        for line in f:
            x, y = line.strip().split()
            # converto in float e salvo in una tupla
            x = float(x)
            y = float(y)
            ghia_data.append((x, y))

        # normalizzo i dati rispetto al massimo valore in modulo
        max_value_ghia = max([abs(y) for x, y in ghia_data])
        ghia_data = [(x, y / max_value_ghia) for x, y in ghia_data]

    # leggo i dati della simulazione dal file
    with open(file_name, "rb") as f:
        # reading velocities
        data = np.fromfile(f, dtype=np.float64)

        max_value = np.abs(data).max()
        if max_value > 0:
            data = data / max_value

        # max_value = max([abs(v) for v in data])
        # data = [v / max_value for v in data]

    return data, ghia_data


def create_graph(data, ghia_data):
    # creo un grafico per i valori di data e lo mostro
    plt.figure()
    # data parte da x=0 a x=1 con passo 1/(len(data)-1)
    N = len(data)
    x = (np.arange(N) + 0.5) / N
    plt.plot(x, data, label="LBM Simulation")

    # Aggiungo i dati di Ghia et al. 1982 per confronto
    ghia_data = np.array(ghia_data)
    plt.plot(ghia_data[:, 0], ghia_data[:, 1], label="Ghia et al.", marker="s")
    plt.xlabel("y coordinate")
    plt.ylabel("u velocity")
    plt.title("Comparison of u-velocity profiles")
    plt.legend()
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <input_file_1> <ghia_bench_data>")
        sys.exit()

    data, ghia_data = read_data(sys.argv[1], sys.argv[2])
    create_graph(data, ghia_data)
