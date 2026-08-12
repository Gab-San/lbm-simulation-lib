import numpy as np
import matplotlib.pyplot as plt

# input file name
input_file = 'vel_norms.txt'


def read_data(file_name):
    with open(file_name, 'r') as f:
        # reading grid size
        nx = int(f.readline().strip())
        ny = int(f.readline().strip())

        # reading velocities
        data = []
        for line in f:
            data.append(float(line.strip()))

    return nx, ny, data

def create_frames(nx, ny, data, num_iterations, vmax):
    frames = []
    #for iter in range(num_iterations):
    iter = 10 # Da 0 a 99

    # legge dati per ogni iterazione e salva in una matrice (ny, nx)
    frame_data = np.array(data[iter * nx * ny:(iter + 1) * nx * ny]).reshape(ny, nx) #.transpose()
    frame_data = frame_data.T   
    
    
    plt.imshow(frame_data, cmap='RdBu_r', origin='lower', vmin=0, vmax = vmax)  # 'origin' è impostato su 'lower' per far partire y da 0 in basso
    plt.colorbar(label='Velocity Magnitude')
    plt.title(f'Iteration {(iter)}')
    plt.show()

    # salva il frame come file immagine
    frame_filename = f'frame_{iter:04d}.png'
    plt.imsave(frame_filename, frame_data, cmap='RdBu_r', vmin=0, vmax=vmax)


    print()  # Per assicurarsi che il prompt successivo inizi su una nuova linea
    return frames

if __name__ == '__main__':
    nx, ny, data = read_data(input_file)
    #print("Stampo la norma delle velocità:\n", data)
    vmax = max(data)
    num_iterations = 100
    frames = create_frames(nx, ny, data, num_iterations, vmax)
