import sys
import matplotlib.pyplot as plt
import numpy as np


def read_ghia_data(bench_name):
    # Leggo i dati del benchmark di Ghia et al.
    ghia_data = []
    with open(bench_name, "r") as f:
        for line in f:
            line_str = line.strip()
            if not line_str:
                continue
            x, y = line_str.split()
            ghia_data.append((float(x), float(y)))

    # Normalizzo i dati di Ghia rispetto al massimo valore in modulo
    max_value_ghia = max([abs(y) for _, y in ghia_data])
    if max_value_ghia > 0:
        ghia_data = [(x, y / max_value_ghia) for x, y in ghia_data]

    return ghia_data


def read_lbm_data(file_name):
    # Leggo i dati della simulazione LBM dal file binario (float64 raw)
    with open(file_name, "rb") as f:
        data = np.fromfile(f, dtype=np.float64)

    # Normalizzo i dati LBM
    max_value = np.abs(data).max() if len(data) > 0 else 0
    if max_value > 0:
        data = data / max_value

    return data


def create_graph(bgk_data, trt_data, ghia_data, output_file=None):
    plt.figure(figsize=(8, 6))

    # Profili LBM da 0 a 1, un colore diverso per ciascun operatore di collisione
    x_bgk = np.linspace(0, 1, len(bgk_data))
    plt.plot(x_bgk, bgk_data, label="BGK", color="tab:blue", linewidth=2)

    x_trt = np.linspace(0, 1, len(trt_data))
    plt.plot(x_trt, trt_data, label="TRT", color="tab:orange", linewidth=2)

    # Dati di confronto Ghia et al. 1982
    ghia_arr = np.array(ghia_data)
    plt.plot(
        ghia_arr[:, 0],
        ghia_arr[:, 1],
        label="Ghia et al. (1982)",
        color="tab:red",
        marker="s",
        linestyle="None",
    )

    plt.xlabel("y coordinate")
    plt.ylabel("u velocity (normalized)")
    plt.title("Comparison of u-velocity profiles")
    plt.legend(fontsize="small", loc="best")
    plt.grid(True)

    # Salva il grafico se viene specificato un percorso di output
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"Grafico salvato in: {output_file}")

    plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(
            f"Usage: {sys.argv[0]} <input_file_bgk.txt> <input_file_trt.txt> "
            f"<ghia_bench_data.txt> [output_image.png]"
        )
        sys.exit(1)

    file_bgk = sys.argv[1]
    file_trt = sys.argv[2]
    file_ghia = sys.argv[3]
    out_img = sys.argv[4] if len(sys.argv) > 4 else "u_velocity_profile.png"

    bgk_data = read_lbm_data(file_bgk)
    trt_data = read_lbm_data(file_trt)
    ghia_data = read_ghia_data(file_ghia)

    create_graph(bgk_data, trt_data, ghia_data, output_file=out_img)
