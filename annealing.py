"""
A module to be used as a template when trying new optimization methods.
"""
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

import earth_data
import geometry
from genetics import random_population
from timer import timeit
import constants
import plot_funcs

dtime = datetime(2015, 5, 1)
clear_polys, big_clear_poly = earth_data.get_clear(dtime)
clear_points = earth_data.get_clear_polys(dtime, just_coords=True)
initial_points = random_population(clear_polys)[0]

T_init = 10       # Initial temperature
T_final = 0.001   # Final temperature
T_steps = 200      # Number of different temperatures in the cooling sequence
T_spin_flips = 100    # Number of spin flips at each temperature
T_schedule = np.linspace(T_init, T_final, T_steps) # Linear cooling schedule
num_anneals = 50  # number of times to run the annealing procedure


def cost(points):
    """
    Calculate the cost of a set of points,
    at a time where the clear area on Earth is represented by big_clear_poly
    """
    return geometry.fov_coverage([points], big_clear_poly)


def energy(points):
    """
    Compute the energy of a set of observations.

    Params:
        :points: A 1D numpy array of (lon, lat) tuples, points to observe
    Return:
        The energy of the configuration if we say cost = energy.
    """
    return 1000 / np.sum(cost(points), axis=1)


def metropolis_decision(initial_energy, flipped_energy, T):
    """
    Use the Metropolis criterion to decide whether to accept a new state. 

    Params:
        :initial_energy: energy of starting configuration
        :flipped_energy: energy of another configuration we are comparing with
        :T: the current temperature of the system
    
    Return:
        True if we should accept the flipped_energy state False otherwise
    """
    if flipped_energy < initial_energy:
        return True
    else:
        acceptance_prob = np.exp((initial_energy - flipped_energy) / T)

        if np.random.rand() < acceptance_prob:
            return True
        else:
            return False


def get_radius(T):
    """get radius of how far away to move a point during SA"""
    return constants.re_km * (T / T_init)


def move_point_within(point, radius):
    # all possible points to move to are given by clear poly centers
    # find all points to move to within the radius
    possible_points = []
    for c in clear_points:
        lon1, lat1 = point
        lon2, lat2 = c
        r = geometry.arc_length(lon1, lat1, lon2=lon2, lat2=lat2)
        if r < radius:
            possible_points.append(c)
    rand_index = int(np.random.randint(0, high=len(possible_points)))
    return possible_points[rand_index]


def metropolis_sa(points, point_to_move, T):
    """
    Run a single Monte Carlo step of SA. Flip a spin and decide whether to accept 
    the new configuration. Call the metropolis_decision method to help you here!
    
    Params:
        :points: A 1D numpy array of (lon, lat) tuples, points to observe.
        :point_to_move: An integer indicating which point should be moved.
        :T: the current temperature of the system

    Return:
        :return: The updated (or unchanged) lattice.
    """
    # Copy the lattice
    points_moved = points.copy()
    
    # Move the point to some other point within a certain radius
    point = points[point_to_move]
    radius = get_radius(T)
    points_moved[point_to_move] = move_point_within(point, radius)
    
    # Decide whether or not to accept
    initial_energy = energy(points)
    flipped_energy = energy(points_moved)
    if metropolis_decision(initial_energy, flipped_energy, T):
        return points_moved
    else:
        return points


def simulated_annealing(T_schedule):
    """
    Run simulated annealing.
    
    Params:
        :T_schedule: A numpy array that is the list of temperatures.

    Return:
        :return: A list containing the energy of the system at the end
                 of each temperature step.
        :return: The final version of the lattice.
    """
    # Keeps track of the final energy at every temperature
    energy_per_step = np.zeros(len(T_schedule))

    # Initialize a random starting point for the lattice
    points = initial_points

    # Work through the temperature schedule and at each point
    #  - Perform T_sweeps attempts to flip random spins on the lattice
    #  - Store the energy of the final configuration in energy_per_step
    for step in range(len(T_schedule)):
        print("Now working on step", step+1, "out of", len(T_schedule))
        T = T_schedule[step]

        for flip in range(T_spin_flips):
            print(f"Flip = {flip:02}\r", end="")
            point_to_move = np.random.randint(0, len(points))
            points = metropolis_sa(points, point_to_move, T)
        plot_funcs.plot_points(points, "a"+str(step), dtime, show=False)
        energy_per_step[step] = energy(points)

    return energy_per_step, points


# Let's plot the energies...
energy_per_step, points = simulated_annealing(T_schedule)
#plt.plot(energy_per_step)
#plt.savefig("e_per_step.png")
plot_funcs.plot_points(points, "Annealing_g0", dtime, show=False)
print(f"Minimum energy value found is {np.min(energy_per_step):.4f}")
best_energies = np.zeros(num_anneals)

for anneal_run in range(num_anneals):
    energy_per_step, points = simulated_annealing(T_schedule)
    plot_funcs.plot_points(points, f"Annealing_g{anneal_run+1}", dtime, show=False)
    best_energies[anneal_run] = energy_per_step[-1]

# Plot a histogram of the energy outputs
# We will store our best energy for comparison later
our_best_energy = np.min(best_energies)
plt.hist(best_energies)
plt.show()
print(f"Minimum energy found by SA is {our_best_energy:.4f}")
