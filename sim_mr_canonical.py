__author__ = 'hok1'

import numpy as np
from multiprocessing import Pool

def raw_sim_particle_levels(N, totalE, R=float('inf')):
    """
    Simulate the energy levels for all particles, given
    there are N particles, totalE unit of energies and R
    energy levels.

    Args:
        N (int): Total number of particles.
        totalE (int): Total units of energies.
        R (int or float, optional): Number of energy levels. Defaults to float('inf').

    Returns:
        ndarray: An array (ndarray) of N elements storing the energy levels of all the
        N particles in the simulation.
    """
    partlevels = np.zeros(N)
    for i in range(totalE):
        avail_part = np.where(partlevels<(R-1))[0]
        #idx = np.random.choice(avail_part)
        idx = avail_part[np.random.randint(0, high=len(avail_part))]
        partlevels[idx] += 1
    return partlevels

def stat_collecting(partlevels):
    """
    Count the number of particles in each energy levels.

    Args:
        partlevels (ndarray): Array of integers indicating the energy levels of each
        particles.

    Returns:
        list of tuples: A list of tuples, where  each tuple includes the energy level
        and the number of particles in that level. For example, it returns

         [(0, 999800), (1, 198), (2, 2)]

        meaning in level 0, there are 999800 particles, level 1, 198 and level 2, 2.
    """
    levels = np.array(sorted(np.unique(partlevels)))
    degeneracies = map(lambda level: len(np.where(partlevels==level)[0]), levels)
    return zip(levels, degeneracies)

def simulate_onethread((N, totalE, R)):
    return stat_collecting(raw_sim_particle_levels(N, totalE, R))

class CanonicalEnsembleSimulationWrapper:
    """
    This is a wrapper class for packaging the simulation using the functions
    @raw_sim_particle_levels and @stat_collecting, while employing parallel
    calculation using the MapReduce method.

    Users have to create an instance of this class, but no arguments are
    needed in the constructor.
    """
    def reduce_simulation(self, level_degeneracy_pairs_lists):
        """
        This reduces the multiple lists of level-degeneracy tuples
        output from a pool of simulation workers into a dictionary.
        For example, if the following two lists are given as follow:
            [(0, 100), (1, 10)] and [(0, 120), (1, 12), (2, 1)],
        it will combine the two lists of tuples and output the dictionary
        that sums up the degeneracies as follow:
            {0: 220, 1: 22, 2: 1}.

        Args:
            level_degeneracy_pairs_lists (list of lists of tuples): Multiple
            lists of tuples, where each tuple contains the energy level and
            degeneracy.

        Returns:
            dict: Dictionary that stores the sums of degeneracies for all
            energy levels.
        """
        occupancies = {}
        for pairs in level_degeneracy_pairs_lists:
            for level, degeneracy in pairs:
                if occupancies.has_key(level):
                    occupancies[level] += degeneracy
                else:
                    occupancies[level] = degeneracy
        return occupancies

    def simulate_multithread(self, N, totalE, R=float('inf'), numthreads=1):
        """
        Simulate, in parallel, the energy levels for all particles, given
        there are N particles, totalE unit of energies and R energy levels.

        Args:
            N (int): Total number of particles.
            totalE (int): Total units of energies.
            R (int or float, optional): Number of energy levels. Defaults to float('inf').
            numthreads (int, optional): Number of workers in the pool. Defaults to 1.

        Returns:
            dict: Dictionary that stores the sums of degeneracies for all
            energy levels.
        """
        if numthreads==1:
            level_degeneracies_pairs_lists = [simulate_onethread((N, totalE, R))]
        else :
            p = Pool(processes=numthreads)
            params = [(N/numthreads, totalE/numthreads, R)]*(numthreads-1)
            params.append((N/numthreads+N % numthreads, totalE/numthreads + totalE % numthreads, R))
            level_degeneracies_pairs_lists = p.map(simulate_onethread, params)
        return self.reduce_simulation(level_degeneracies_pairs_lists)

    def stat_analyze(self, occupancies):
        """
        Analyze the energy levels and their degeneracies. Outputs
        the fitting coefficients in the exponential distribution
        alpha and -beta, where the distribution ~ exp(alpha-beta*n),
        where n is the energy level, the average energy and its error
        (standard deviation).

        Args:
            occupancies (dict): Dictionary that stores the sums of degeneracies for all
            energy levels.

        Returns:
            ndarray: Array of floats with two elements. The first is -beta, and
            the second is alpha.
            float: Average energy.
            float: Error (standard deviation) of energy.
        """
        levels = sorted(occupancies.keys())
        degeneracies = map(lambda level: occupancies[level], levels)

        levels = np.array(levels)
        degeneracies = np.array(degeneracies)
        coefs = np.polyfit(levels, np.log(degeneracies), 1)

        mean_energy = np.sum(levels*degeneracies) / np.sum(degeneracies)
        std_energy = np.sqrt(np.sum((levels-mean_energy)*(levels-mean_energy)*degeneracies) / np.sum(degeneracies))
        return coefs, mean_energy, std_energy

def main():
    """
    Main block of the program.
    """
    N = int(raw_input("Number of particles? "))
    R = float(raw_input("Number of energy levels? "))
    totalE = int(raw_input("Total energy? "))
    num_workers = int(raw_input("Number of workers? "))

    simulator = CanonicalEnsembleSimulationWrapper()

    occupancies = simulator.simulate_multithread(N, totalE, R, numthreads=num_workers)

    print 'Results:'
    print 'level\tdegeneracy'
    for level in sorted(occupancies.keys()):
        print level, '\t', occupancies[level]

    coefs, mean_energy, std_energy = simulator.stat_analyze(occupancies)
    print 'alpha = ', coefs[1], ', beta = ', -coefs[0]
    print 'average energy = ', mean_energy
    print 'standard deviation = ', std_energy

if __name__ == '__main__':
    main()