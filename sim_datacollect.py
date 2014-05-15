__author__ = 'hok1'

import sim_mr_canonical as sim
import csv

N = 100000000
total_energies = [10, 50, 100, 500, 1000, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000]
num_workers = 10
nrep = 20
filename = 'raw_canonical_sim.csv'

if __name__ == '__main__':
    simulator = sim.CanonicalEnsembleSimulationWrapper()
    outf = open(filename, 'wb')
    writer = csv.writer(outf)
    header = ['N', 'total_energy', 'mean_energy', 'std_energy', 'beta', 'alpha']
    writer.writerow(header)
    for total_energy in total_energies:
        for i in range(nrep):
            occupancies = simulator.simulate_multithread(N, total_energy, numthreads=num_workers)
            coefs, mean_energy, std_energy = simulator.stat_analyze(occupancies)
            writer.writerow([N, total_energy, mean_energy, std_energy, -coefs[0], coefs[1]])
    outf.close()