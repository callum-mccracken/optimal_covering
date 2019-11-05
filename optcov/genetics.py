import numpy as np
from shapely.geometry import Polygon, MultiPolygon, Point
from shapely.errors import TopologicalError
import itertools
import multiprocessing

from optcov.utils.timer import timeit
from optcov.utils import geometry
from optcov.utils import sampling
from optcov import plot_funcs, earth_data

num_generations = 10
num_parents = 10
num_members = 20
num_fovs = 30
num_mutations = 5

pop_size = (num_members, num_fovs)


@timeit
def random_population(clear_polys):
    """will generate random fovs, tuples of the form (lon, lat)
    (but not super random, they should at least be clear spots!)
    """
    # we'll return this at the end
    fovs = np.empty(pop_size, dtype=tuple)

    # some points where the sky is clear, visible, etc.
    # note clear_coords are in (lon, lat) format
    clear_coords = []
    for p in clear_polys:
        lon, lat = p.centroid.xy  # returns these as odd arrays for some reason
        lon, lat = float(lon[0]), float(lat[0])
        clear_coords.append((lon, lat))

    # shuffle them up since we're going to pick a few at random
    np.random.shuffle(clear_coords)

    # pick the first few random ones
    # - I assume we have way more clear spots than fovs, so this is safe
    for index, _ in np.ndenumerate(fovs):
        coords = clear_coords.pop()
        fovs[index] = coords
    return fovs


@timeit
def calc_fitness(population, clear_polys, fitness_type="coverage"):
    """General fitness function to be maximized.

    This function returns an array containing the fitness of each member
    of the population.

    - population is a 2d array of (lon, lat) tuples.

    - fitness_type is a string --> which type of fitness calculation to do

    - mask_polys is a variable required to calculate fitness
      for the "coverage" fitness_type.

    - feel free to add your own elif cases
      if you want to add extra fitness types
    """
    if fitness_type == "coverage":
        fov_fitness = geometry.fov_coverage(population, clear_polys)
    else:
        raise ValueError("invalid fitness_type" + str(fitness_type))
    # this returns the fitness at the level of fovs, if you're looking
    # for the fitness of members, call this function then use this line after:
    # member_fitness = np.sum(fov_fitness, axis=1)
    return fov_fitness


@timeit
def select_parents(population, member_fitness):
    # Select the best individuals in the current generation as parents
    # for producing the offspring of the next generation.

    # reorder the population in terms of fitness, largest to smallest
    fitness_order = list(reversed(member_fitness.argsort()))
    parents = [population[f, :] for f in fitness_order[:num_parents]]
    parental_fitness = [member_fitness[f] for f in fitness_order[:num_parents]]
    return parents, parental_fitness


@timeit
def offspring(parents, parental_fitness, clear_polys,
              algorithm="best_combination"):
    offspring_size = (num_members-num_parents, num_fovs)
    offspring = np.empty(offspring_size, dtype=tuple)

    if algorithm == 'crossover':
        # The point at which crossover takes place between two parents.
        # Usually, it is at the center.
        crossover_point = np.uint8(offspring_size[1]/2)

        for i in range(len(offspring)):
            # Index of the first parent to mate.
            parent1_idx = i % len(parents)
            # Index of the second parent to mate.
            parent2_idx = (i+1) % len(parents)

            if parental_fitness is None:
                # The new offspring will have its first few
                # fovs taken from the first parent.
                offspring[i, 0:crossover_point] = parents[
                    parent1_idx, 0:crossover_point]
                # The rest will be taken from the second parent.
                offspring[i, crossover_point:] = parents[
                    parent2_idx, crossover_point:]
            else:
                # same idea as above, but we'll take the best fovs from each

                # find best fovs in each parent:
                parent1 = parents[parent1_idx]
                parent2 = parents[parent2_idx]
                p1_fitness = parental_fitness[parent1_idx, :]
                p2_fitness = parental_fitness[parent2_idx, :]
                # sort both arrays in the same way, i.e. the order of fitness
                # https://stackoverflow.com/questions/1903462/how-can-i-zip-sort-parallel-numpy-arrays
                p1_order = p1_fitness.argsort()
                p2_order = p2_fitness.argsort()
                parent1 = parent1[p1_order]
                parent2 = parent2[p2_order]
                # reverse parent1 so when you take the first few you get
                # the largest fitness, not smallest
                parent1 = np.flip(parent1)
                offspring[i, 0:crossover_point] = parent1[0:crossover_point]
                # The rest will be taken from the second parent.
                offspring[i, crossover_point:] = parent2[crossover_point:]
    elif algorithm == 'best_combination':
        # take possible combinations of parents
        # num_fovs = offspring_size[1]
        # all possible combinations of parents (may be large, not HUGE)
        parent_pairings = list(itertools.combinations(parents, 2))
        fov_combos = []
        # there will be (n_parents choose 2) pairings, with n samples
        # we'll take n combinations per pairing, and grab the best one later
        n = 10
        for pairing in parent_pairings:
            parent1, parent2 = pairing
            # take n of the possible arrangements of the fovs from parents
            # if we took all possible arrangemland_maskat would be
            # 2n_fovs choose n_fovs, which is not feasible
            fov_combos.append(sampling.sampled(parent1, parent2, n))

        fov_combos = np.array(fov_combos)

        print("checking", len(fov_combos), "pairings")
        print("each with", n, "arrangements")
        # this thing is a calculation of population fitness for all our
        # n_parents choose 2 pairings
        combos_fitness = []
        for c in fov_combos:
            try:
                # sometimes this just doesn't work, not sure why...
                combo_fitness = calc_fitness(c, earth_data.big_clear_poly)
            except TopologicalError:
                combo_fitness = 0
            combos_fitness.append(combo_fitness)


        # now we must pick the best arrangement (or combo) in each pairing
        best_combos = []
        fitnesses = []
        for pair_index, pair_fitnesses in enumerate(combos_fitness):
            fovs = parent_pairings[pair_index]
            # do a sum, so we get the total fitness of each combo,
            # not just the fitness of each individual FoV
            sum_fitnesses = np.sum(pair_fitnesses, axis=1)
            # sort/flip so theclear_maskndex is greatest fitness, last is least
            fitness_order = sum_fitnesses.argsort()
            sum_fitnesses = np.flip(sum_fitnesses[fitness_order])
            fovs = np.flip(fovs[fitness_order])
            for i in range(len(sum_fitnesses)):
                best_combos.append(fovs)
                fitnesses.append(sum_fitnesses)
        best_combos = np.array(best_combos)
        fitnesses = np.array(fitnesses)
        # now best_combos[i][j] is the j'th best arrangement
        # of the i'th pairing, and fitnesses[i][j] is corresponding fitness
        # so now set the offspring
        for i, _ in enumerate(offspring):
            offspring[i] = best_combos[i][0]  # best arrangement of ith pair
    elif algorithm == 'just return':
        return parents
    return offspring


@timeit
def mutation(offspring, member_fitness, clear_polys):
    # Mutation changes members in each offspring collection randomly.
    for i, member in enumerate(offspring):
        fov_fitness = member_fitness[i]
        fovs = offspring[i]
        # get order of fitness, largest to smallest
        fov_fitness_order = reversed(list(fov_fitness.argsort()))
        # re-order offspring and fovs in terms of fitness order
        obs_polying_fitness = member_fitness[fov_fitness_order]
        fovs = fovs[fov_fitness_order]
        # get clear points on Earth
        clear_coords = [p.centroid for p in clear_polys]
        # randomize since we'll use these for random mutations
        np.random.shuffle(clear_coords)

        # now mutate the worst 'num_mutations' fovs
        for fov_number in range(num_mutations):
            fov_to_change = member[-1*fov_number]  # take from the end of list

            other_fovs = [f for f in fovs if f is not fov_to_change]

            assert type(fov_to_change) == Polygon

            centroid = fov_to_change.centroid

            # now we're going to pick a new spot to put the fov.

            # ensure this new point is not inside of any other FoVs
            # some overlap will still be possible, but the center of
            # one fov should logically never be inside another fov

            in_others = True
            while in_others:
                # pick a random clear spot
                new_lon, new_lat = clear_coords.pop()

                # assume not in any others, then check
                in_others = False
                for other_f in other_fovs:
                    if other_f.contains(Point(new_lon, new_lat)):
                        in_others = True
                        break

                # if we run out of points to pick, we'll raise an error
                # when calling pop(), but that's never been an issue before

            # make polygon from that new point
            vertices = geometry.generate_vertices((new_lon, new_lat))
            mutated_fov = Polygon(vertices)

            # ensure it's valid
            if not mutated_fov.is_valid:
                raise ValueError("invalid")

            # set the fov of the offspring to the new mutated fov
            member[fov_number] = mutated_fov
    return offspring


def perform_genetics(dtime):
    """
    This is the function which does the main process in the genetic algorithm
    - generate initial population
    - evolve
    - repeat until you're done a number of iterations
    """
    # get clear polygons at the given time
    print("getting clear polys")
    clear_polys = earth_data.clear_polys
    clear_poly = earth_data.big_clear_poly
    print("plotting clear area")
    # plot_funcs.plot_clear(dtime, clear_polys=clear_polys)

    # generate population
    print("generating population")
    population = random_population(clear_polys)

    best_fitnesses = []  # to store best fitness of each generation

    # now iterate and 'evolve'
    for generation in range(num_generations):
        print("generation", generation)
        # calculate fitness
        pop_fitness = calc_fitness(population, clear_poly)
        member_fitness = np.sum(pop_fitness, axis=1)

        # Find the current best member
        best_fitness = np.max(member_fitness)
        best_fitnesses.append(best_fitness)
        # if multiple maxes, take first option
        best_match_idx = list(member_fitness).index(best_fitness)
        best_member = population[best_match_idx]
        print("Best fitness:", best_fitness, "index", best_match_idx)

        # plot the best member
        # plot_funcs.plot_member(best_member, generation, clear_polys, dtime)
        parents, parental_fitness = select_parents(population, member_fitness)

        # Generate the next generation using parents
        kids = offspring(parents, parental_fitness, clear_polys)

        # Add some variations to the offspring, i.e. mutations
        offspring_mutation = mutation(kids, pop_fitness, clear_polys)

        # Create the new population based on the parents and offspring
        population[0:num_parents] = parents
        population[num_parents:] = offspring_mutation

    # Get the best solution after finishing all generations.
    print("Generation : ", generation+1)
    pop_fitness = calc_fitness(population, clear_polys)
    member_fitness = np.sum(pop_fitness, axis=1)
    best_fitness = np.max(member_fitness)
    # if multiple maxes, take first option
    best_match_idx = list(member_fitness).index(best_fitness)
    best_member = population[best_match_idx]
    best_fitnesses.append(best_fitness)

    return population[best_match_idx]  # <- array of 'best' FOVs

if __name__ == "__main__":
    import optcov.utils.constants as c
    perform_genetics(c.dtime)
