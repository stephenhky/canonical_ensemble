__author__ = 'hok1'

import numpy as np
from multiprocessing import Pool

def raw_sim_particle_levels((N, R, totalE)):
    partlevels = np.zeros(N)
    for i in range(totalE):
        avail_part = np.where(partlevels<(R-1))[0]
        idx = np.random.choice(avail_part)
        partlevels[idx] += 1
    return partlevels

def sim_particle_levels(N, R, totalE, num_workers=1):
    p = Pool(processes=num_workers)
    params = [(N/num_workers, R, totalE/num_workers)]*num_workers
    part_levels = p.map(raw_sim_particle_levels, params)
    return np.concatenate(part_levels)

def stat_collecting(partlevels):
    levels = np.array(sorted(np.unique(partlevels)))
    degeneracies = map(lambda level: len(np.where(partlevels==level)[0]), levels)
    return levels, degeneracies

def stat_analyze(partlevels, toprint=False):
    levels, degeneracies = stat_collecting(partlevels)
    coefs = np.polyfit(levels, np.log(degeneracies), 1)
    mean_energy = np.mean(partlevels)
    std_energy = np.std(partlevels)
    if toprint:
        print 'alpha = ', coefs[1], ', beta = ', coefs[0]
        print 'average energy = ', mean_energy
        print 'standard deviation = ', std_energy
    return coefs, mean_energy, std_energy

def main():
    N = int(raw_input("Number of particles? "))
    R = int(raw_input("Number of energy levels? "))
    totalE = int(raw_input("Total energy? "))
    num_workers = int(raw_input("Number of workers? "))
    stat_analyze(sim_particle_levels(N, R, totalE, num_workers=num_workers), toprint=True)

if __name__ == '__main__':
    main()