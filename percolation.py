import numpy as np
from scipy.special import expit
import random
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Define the size of the lattice and the number of iterations
N = 50
iterations = 1000

# Define the parameters of the simulation
J = 1
H = -10
kT = 0.01

# Define the lattice as an array of random spins (+1 or -1)
lattice = np.random.choice([-1, 1], size=(N, N))

# Define a function to calculate the energy of a site
def energy(site, JJ, HH):
    m_i = np.sum(lattice[site[0], site[1]] * np.array([
        lattice[(site[0]-1)%N, site[1]],
        lattice[(site[0]+1)%N, site[1]],
        lattice[site[0], (site[1]-1)%N],
        lattice[site[0], (site[1]+1)%N]
    ]))
    return -JJ*m_i - HH, +JJ*m_i + HH

# Define a function to flip a spin according to the Metropolis algorithm
def flip_spin(site):
    E_plus, E_minus = energy(site, J, H)
    prob = expit(-(E_minus - E_plus) / kT)
    # print(prob)
    if random.uniform(0, 1) < prob:
        lattice[site[0], site[1]] = 1
    else:
        lattice[site[0], site[1]] = -1

# Define a function to update the plot
def update(frame):
    for i in range(iterations):
        site = (random.randint(0, N-1), random.randint(0, N-1))
        flip_spin(site)
    im.set_array(lattice)
    return [im]

if __name__ == "__main__":

    # Create the plot and animation
    fig, ax = plt.subplots()
    im = ax.imshow(lattice, cmap='binary')
    ani = animation.FuncAnimation(fig, update, frames=100, interval=50, blit=True)

    plt.show()
