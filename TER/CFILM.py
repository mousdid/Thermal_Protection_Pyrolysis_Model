import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.animation import FuncAnimation

def read_temperature_data(filename):
    try:
        with open(filename, 'r') as sauv:
            return [float(line) for line in sauv.readlines()]
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return None

def update(frame):
    global Tgrid, colormap, norm, map_ax, colorbar

    temperature_filename = f'{frame}_RHO.txt'
    T = read_temperature_data(temperature_filename)
    
    if T is not None:
        Tgrid = np.reshape(T, (Nx, Ny))

        # Update the image data
        img.set_array(Tgrid)

        # Set clim based on the new data
        img.set_clim(vmin=Tgrid.min(), vmax=Tgrid.max())
        
        # Update colorbar label
        colorbar.set_label(f'Echelle(K) - Frame {frame}')
    else:
        print(f"Error: Data not available for frame {frame}")

with open('data.txt', 'r') as file:
    Nx = int(file.readline().strip().split("=")[1])
    Ny = int(file.readline().strip().split("=")[1])

initial_frame = 0
Tgrid = np.reshape(read_temperature_data(f'{initial_frame}_RHO.txt'), (Nx, Ny))
colormap = plt.cm.jet
norm = Normalize(vmin=Tgrid.min(), vmax=Tgrid.max())

fig1, (map_ax) = plt.subplots(figsize=(8, 6))

img = map_ax.imshow(Tgrid, cmap=colormap, norm=norm, interpolation='bilinear', extent=[0, Nx, 0, Ny])
colorbar = fig1.colorbar(img, ax=map_ax, label='Echelle(kg/m3)')
map_ax.set_title('Masse volumique')
map_ax.set_xlabel('X')
map_ax.set_ylabel('Y')

# Disable cursor interaction for the imshow plot
#map_ax.set_cursor_props(visible=False)

animation_frames = range(initial_frame, 100, 1)
animation = FuncAnimation(fig1, update, frames=animation_frames, interval=10, repeat=True)

with open('temp_mailles', 'r') as mailles:
    maille_data = np.genfromtxt(mailles)  # Assuming data is space-separated

fig2, curves_ax = plt.subplots(figsize=(8, 6))

# Create a new figure and axis for the three curves
curves_ax.set_title('Température en fonction du temps')
curves_ax.set_xlabel('Temps')
curves_ax.set_ylabel('Température')

# Plot the three curves on the same axis
for i in range(maille_data.shape[1] // 3):  # Assuming 3 curves per line
    y_values = maille_data[:, i * 3:i * 3 + 3]
    x_values = np.arange(maille_data.shape[0])  # Use the index as x-values
    curves_ax.plot(x_values, y_values)

legend_labels = [f'Sensor 1',f'Sensor 2',f'Sensor 3']
curves_ax.legend(labels=legend_labels)

plt.show()
