__author__ = 'hok1'

import numpy as np
from multiprocessing import Pool

def raw_sim_particle_levels(N, totalE, R=float('inf')):
    partlevels = np.zeros(N)
    for i in range(totalE):
        avail_part = np.where(partlevels<(R-1))[0]
        idx = np.random.choice(avail_part)
        partlevels[idx] += 1
    return partlevels

def stat_collecting(partlevels):
    levels = np.array(sorted(np.unique(partlevels)))
    degeneracies = map(lambda level: len(np.where(partlevels==level)[0]), levels)
    return zip(levels, degeneracies)

def simulate_onethread((N, totalE, R)):
    return stat_collecting(raw_sim_particle_levels(N, totalE, R))

class CanonicalEnsembleSimulationWrapper:
    def reduce_simulation(self, level_degeneracy_pairs_lists):
        occupancies = {}
        for pairs in level_degeneracy_pairs_lists:
            for level, degeneracy in pairs:
                if occupancies.has_key(level):
                    occupancies[level] += degeneracy
                else:
                    occupancies[level] = degeneracy
        return occupancies

    def simulate_multithread(self, N, totalE, R=float('inf'), numthreads=1):
        if numthreads==1:
            level_degeneracies_pairs_lists = [simulate_onethread((N, totalE, R))]
        else :
            p = Pool(processes=numthreads)
            params = [(N/numthreads, totalE/numthreads, R)]*(numthreads-1)
            params.append((N/numthreads+N % numthreads, totalE/numthreads + totalE % numthreads, R))
            level_degeneracies_pairs_lists = p.map(simulate_onethread, params)
        return self.reduce_simulation(level_degeneracies_pairs_lists)

    def stat_analyze(self, occupancies):
        levels = sorted(occupancies.keys())
        degeneracies = map(lambda level: occupancies[level], levels)

        levels = np.array(levels)
        degeneracies = np.array(degeneracies)
        coefs = np.polyfit(levels, np.log(degeneracies), 1)

        mean_energy = np.sum(levels*degeneracies) / np.sum(degeneracies)
        std_energy = np.sqrt(np.sum((levels-mean_energy)*(levels-mean_energy)*degeneracies) / np.sum(degeneracies))
        return coefs, mean_energy, std_energy

def main():
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
    print 'alpha = ', coefs[1], ', beta = ', coefs[0]
    print 'average energy = ', mean_energy
    print 'standard deviation = ', std_energy

if __name__ == '__main__':
    main()