import numpy as np
import random
from numba import jit, prange
from lib import algo

# @jit(nopython=True)
def fitness(Garr, Parr, phif, psif, CB, CG, ND2, C1, C2):
    Garr = algo.rr(phif, psif, CB, CG, ND2, C1, C2, Garr, Parr)
    return algo.steric_fast(Garr, Parr)

# @jit(nopython=True)
def create_individual(phisd, psisd):
    print(phisd, psisd)
    phisd = list(phisd)  # Assuming phisd is also a DataFrame
    psisd = list(psisd)
    phi = np.random.uniform(phisd[0], phisd[1])
    psi = np.random.uniform(psisd[0], psisd[1])
    return phi, psi

# @jit(nopython=True)
def create_population(phisd, psisd, population_size):
    population = np.empty((population_size, 2))
    for i in prange(population_size):
        phi, psi = create_individual(phisd, psisd)
        population[i, 0] = phi
        population[i, 1] = psi
    return population

# @jit(nopython=True)
def crossover(parent1, parent2):
    crossover_point = np.random.randint(0, 2)
    child1 = np.empty(2)
    child2 = np.empty(2)
    child1[:crossover_point] = parent1[:crossover_point]
    child1[crossover_point:] = parent2[crossover_point:]
    child2[:crossover_point] = parent2[:crossover_point]
    child2[crossover_point:] = parent1[crossover_point:]
    return child1, child2

# @jit(nopython=True)
def mutate(individual, mutation_rate, phisd, psisd):
    mutated_individual = np.empty(2)
    for i in prange(2):
        if np.random.random() < mutation_rate:
            if i == 0:
                mutated_individual[i] = np.random.uniform(phisd[0], phisd[1])
            else:
                mutated_individual[i] = np.random.uniform(psisd[0], psisd[1])
        else:
            mutated_individual[i] = individual[i]
    return mutated_individual

# @jit(nopython=True)
def evolve_population(population, fitness_fn, CB, CG, ND2, C1, C2, Garr, Parr, phisd, psisd, mutation_rate):
    population_size = len(population)
    new_population = np.empty((population_size, 2))
    ranked_population_indices = np.argsort(np.array([fitness_fn(Garr, Parr, ind[0], ind[1], CB, CG, ND2, C1, C2) for ind in population]))

    for i in prange(population_size // 2):
        parent1_index = ranked_population_indices[np.random.randint(0, population_size // 2)]
        parent2_index = ranked_population_indices[np.random.randint(0, population_size // 2)]
        parent1 = population[parent1_index]
        parent2 = population[parent2_index]
        child1, child2 = crossover(parent1, parent2)
        new_population[2*i] = mutate(child1, mutation_rate, phisd, psisd)
        new_population[2*i+1] = mutate(child2, mutation_rate, phisd, psisd)

    return new_population

# @jit(nopython=True)
def genetic_algorithm_opt(G,CB, CG, ND2, C1, C2, Garr, Parr, phisd, psisd, population_size=100, generations=100, mutation_rate=0.1):
    population = create_population(phisd, psisd, population_size)
    
    for _ in range(generations):
        population = evolve_population(population, fitness, CB, CG, ND2, C1, C2, Garr, Parr, phisd, psisd, mutation_rate)
        algo.recorder(Garr,G,'genetic opt')
    best_individual_index = np.argmin(np.array([fitness(Garr, Parr, ind[0], ind[1], CB, CG, ND2, C1, C2) for ind in population]))
    best_individual = population[best_individual_index]
    phif, psif = best_individual
    r = fitness(Garr, Parr, phif, psif, CB, CG, ND2, C1, C2)

    return phif, psif, r
