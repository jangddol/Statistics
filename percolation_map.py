import numpy as np
from scipy.special import expit
import random
import math
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from percolation import energy
from multiprocessing import Pool, set_start_method

# Set the start method to "spawn"
try:
    set_start_method('spawn')
except RuntimeError:
    pass

# Define the size of the lattice and the number of iterations
N = 50
iterations = 30000

# Define the range of kT and H values to test
kT_range = np.linspace(0.01, 0.15, 100)
H_range = np.linspace(-0.2, 0.2, 100)

# Define the lattice as an array of random spins (+1 or -1)
lattice = np.random.choice([-1, 1], size=(N, N))

# Define a colormap for the final results
# cmap = ListedColormap(['black', 'white'])
cmap = 'viridis'

# Define a function to run the simulation for a single value of kT and H
def simulate(kT, H):
    # Initialize the lattice
    lattice = np.random.choice([-1, 1], size=(N, N))
    
    # Run the simulation until it reaches a steady state
    for n in range(iterations):
        site = (random.randint(0, N-1), random.randint(0, N-1))
        E_plus, E_minus = energy(site, 1, H)
        prob = expit((E_plus - E_minus) / kT)
        if random.uniform(0, 1) < prob:
            lattice[site[0], site[1]] = 1
        else:
            lattice[site[0], site[1]] = -1
    
    # Calculate the ratio of black to white spins
    ratio = np.count_nonzero(lattice == 1) / (N*N)
    return ratio

if __name__ == '__main__':
    # Create a pool of worker processes
    with Pool() as p:
        # Loop over each value of kT and H in parallel
        results = p.starmap(simulate, [(kT, H) for kT in kT_range for H in H_range])
    
    testArray = [(kT, H) for kT in kT_range for H in H_range]
    # print(testArray)
    
    # Reshape the results into a 2D array
    ratios = np.array(results).reshape(len(kT_range), len(H_range))
    ratios = np.rot90(ratios)

    # Set the figure size and plot the final ratios as a 2D map
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    image = ax.imshow(ratios, cmap=cmap, extent=(min(kT_range), max(kT_range), min(H_range), max(H_range)), origin='lower', interpolation='nearest', aspect='auto')
    ax.set_xlabel('kT')
    ax.set_ylabel('H')
    cbar = plt.colorbar(image, ax=ax)
    cbar.set_label('Ratio of White Spins')
    plt.show()
