__author__ = 'hok1'

import numpy as np

def sim_particle_levels(N, R, totalE):
    partlevels = np.zeros(N)
    for i in range(totalE):
        avail_part = np.where(partlevels<(R-1))[0]
        idx = np.random.choice(avail_part)
        partlevels[idx] += 1
    return partlevels

def stat_collecting(partlevels):
    levels = np.array(sorted(np.unique(partlevels)))
    degeneracies = map(lambda level: len(np.where(partlevels==level)[0]), levels)
    return levels, degeneracies

def stat_analyze(partlevels):
    levels, degeneracies = stat_collecting(partlevels)
    coefs = np.polyfit(levels, np.log(degeneracies), 1)
    print 'alpha = ', coefs[1], ', beta = ', coefs[0]
    return coefs