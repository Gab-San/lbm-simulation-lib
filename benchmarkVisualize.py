import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# input file name
input_file_1 = 'vel_y_testing.txt'  # Dati presi dalla simulazione
input_file_2 = 'bench_data.txt' # Dati presi da Ghia et al. 1982

def read_data(file_name):

    # leggo i dati di ghia
    with open(input_file_2, 'r') as f:
        ghia_data = []
        for line in f:
            x, y = line.strip().split()
            # converto in float e salvo in una tupla
            x = float(x)
            y = 0.45 * float(y)
            ghia_data.append((x, y))

            #ghia_data.append(float(line.strip()))

    # leggo i dati della simulazione dal file
    with open(file_name, 'r') as f:
      # reading velocities
        data = []
        for line in f:
            data.append(float(line.strip()))

    return data, ghia_data

def create_graph(data, ghia_data):
    # creo un grafico per i valori di data e lo mostro
    plt.figure()
    # data parte da x=0 a x=1 con passo 1/(len(data)-1)
    x = np.linspace(0, 1, len(data))
    plt.plot(x, data, label='LBM Simulation')

    # Aggiungo i dati di Ghia et al. 1982 per confronto
    ghia_data = np.array(ghia_data)
    plt.plot(ghia_data[:, 0], ghia_data[:, 1], label='Ghia et al.', marker='s')
    plt.xlabel('y coordinate')
    plt.ylabel('u velocity')
    plt.title('Comparison of u-velocity profiles')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    data, ghia_data = read_data(input_file_1)
    create_graph(data, ghia_data)