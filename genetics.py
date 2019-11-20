import numpy as np
from shapely.geometry import Point
import itertools

from timer import timeit
import geometry
import sampling
import plot_funcs

num_generations = 2
num_parents = 10
num_members = 20
num_fovs = 30
num_mutations = 5

pop_size = (num_members, num_fovs)


def random_population(clear_polys):
    """
    this will generate random fovs, tuples of the form (lon, lat)
    (but not super random, they should at least be clear spots!)

    clear_polys is a list of Polygons which define the clear area
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
def calc_fitness(population, big_clear_poly, fitness_type="coverage"):
    """
    General fitness function to be maximized.

    This function returns an array containing the fitness of each member
    of the population.

    - population is a 2d array of (lon, lat) tuples.

    - big_clear_poly is the union of all clear spots on the earth

    - fitness_type is a string --> which type of fitness calculation to do

    - feel free to add your own elif cases if you want extra fitness types
    """
    if fitness_type == "coverage":
        fov_fitness = geometry.fov_coverage(population, big_clear_poly)
    else:
        raise ValueError("invalid fitness_type" + str(fitness_type))
    # this returns the fitness at the level of fovs, if you're looking
    # for the fitness of members, call this function then use this line after:
    # member_fitness = np.sum(fov_fitness, axis=1)
    return fov_fitness


def select_parents(population, member_fitness):
    """Select the best individuals in the current generation as parents
    for producing the offspring of the next generation.
    
    population = list of lists of (lon, lat) tuples

    member_fitness = list of floats, fitnesses of each member of population
    """
    # reorder the population in terms of fitness, largest to smallest
    fitness_order = list(reversed(member_fitness.argsort()))
    # take the num_parents best ones
    parents = [population[f, :] for f in fitness_order[:num_parents]]
    parental_fitness = [member_fitness[f] for f in fitness_order[:num_parents]]
    return parents, parental_fitness


@timeit
def offspring(parents, parental_fitness, big_clear_poly,
              algorithm="best_combination"):
    """
    Get offspring given parents and their fitness

    feel free to add other offpring algorithms, etc.

    parents = slice of population, list of lists of (lon, lat) tuples

    parental_fitness = 2d array of floats, fitness of each parent's points
    """
    offspring_size = (num_members-num_parents, num_fovs)
    offspring_arr = np.empty(offspring_size, dtype=tuple)

    if algorithm == 'best_combination':
        # take possible combinations of parents
        # n of the possible arrangements of the fovs from parents
        # if we took all possible arrangements, we would have
        # 2n_fovs choose n_fovs, which is not feasible to calculate
        n = 10
        # all possible pairs of parents
        pairs = list(itertools.combinations(parents, 2))
        # n possible new members taking some fovs from each parent
        fov_combos = np.array([sampling.sampled(*p, n) for p in pairs])

        print("checking", len(fov_combos), "pairings")
        print("each with", n, "arrangements")

        # calculate fitness of each combination, this takes quite a while
        combos_fitness = [calc_fitness(combo, big_clear_poly)
                          for combo in fov_combos]

        # figure out which combinations were the best
        best_combos = []
        fitnesses = []
        for pair_index, pair_fitnesses in enumerate(combos_fitness):
            # members we generated for this pairing
            members = fov_combos[pair_index]
            # fitnesses of these members
            member_fitnesses = np.sum(pair_fitnesses, axis=1)
            # order of fitness elements from greatest to least
            fitness_order = np.flip(member_fitnesses.argsort())
            # append fovs and fitnesses in fitness order
            best_combos.append([members[i] for i in fitness_order])
            fitnesses.append([member_fitnesses[i] for i in fitness_order])
        best_combos = np.array(best_combos)
        fitnesses = np.array(fitnesses)

        # now best_combos[i][j] is the j'th best arrangement
        # of the i'th pairing, and fitnesses[i][j] is corresponding fitness
        # so now set the offspring
        for i, _ in enumerate(offspring_arr):
            offspring_arr[i] = best_combos[i][0]  # best arrangement of ith pair
    return offspring_arr


@timeit
def mutation(offspring, offspring_fitness, clear_polys):
    """
    Mutation changes members in each offspring collection randomly.
    Note that offspring fovs are already ordered from best to worst!
    """
    for i, member in enumerate(offspring):
        # get clear points on Earth
        clear_coords = [p.centroid for p in clear_polys]
        # randomize since we'll use these for random mutations
        np.random.shuffle(clear_coords)

        # now mutate the worst 'num_mutations' fovs
        for fov_number in range(num_mutations):
            mutate_index = -1*fov_number
            fov_to_change = member[mutate_index]  # take from the end of list

            # the other fovs in the member
            other_fovs = [f for f in member if f is not fov_to_change]

            # now we're going to pick a new spot to put the fov.

            # ensure this new point is not inside of any other FoVs.
            # some overlap will still be possible, but the center of
            # one fov should logically never be inside another fov

            in_others = True
            while in_others:
                # pick a new random clear spot
                clear_point = clear_coords.pop()
                new_lon, new_lat = clear_point.x, clear_point.y

                # assume not in any others, then check
                in_others = False
                for other_f in other_fovs:
                    other_f_poly = geometry.obs_poly(*other_f)
                    if other_f_poly.contains(Point(new_lon, new_lat)):
                        in_others = True
                        break

                # if we run out of points to pick, we'll raise an error
                # when calling pop(), but that's never been an issue before

            # put that mutated fov back in the offspring
            offspring[i][mutate_index] = (new_lon, new_lat)
    return offspring


def do_genetics(dtime, clear_polys, big_clear_poly):
    """
    This is the function which does the main process in the genetic algorithm
    - generate initial population
    - evolve
    - repeat until you're done a number of iterations
    """
    # generate population
    print("generating population")
    population = random_population(clear_polys)

    best_fitnesses = []  # to store best fitness of each generation

    # now iterate and 'evolve'
    for generation in range(num_generations):
        print("generation", generation)

        # calculate fitness, 2d array of floats
        pop_fitness = calc_fitness(population, big_clear_poly)
        member_fitness = np.sum(pop_fitness, axis=1)

        # Find the current best member
        best_fitness = np.max(member_fitness)
        best_fitnesses.append(best_fitness)
        # (if multiple maxes, take first option)
        best_match_idx = list(member_fitness).index(best_fitness)
        best_member = population[best_match_idx]
        print("Best fitness:", best_fitness, "index", best_match_idx)

        # plot the best member
        title = "gen_{}".format(generation)
        plot_funcs.plot_points(best_member, title, dtime, show=False)

        # get parents
        parents, parental_fitness = select_parents(population, member_fitness)

        # generate the next generation using parents
        kids = offspring(parents, parental_fitness, big_clear_poly)

        # mutate the offspring
        offspring_mutation = mutation(kids, pop_fitness, clear_polys)

        # create the new population based on the parents and offspring
        population[0:num_parents] = parents
        population[num_parents:] = offspring_mutation

    # get the best solution after finishing all generations.
    print("Generation : ", generation+1)
    pop_fitness = calc_fitness(population, big_clear_poly)
    member_fitness = np.sum(pop_fitness, axis=1)
    best_fitness = np.max(member_fitness)
    # (if multiple maxes, take first option)
    best_match_idx = list(member_fitness).index(best_fitness)
    best_member = population[best_match_idx]
    best_fitnesses.append(best_fitness)

    return population[best_match_idx]  # <- array of 'best' points to look