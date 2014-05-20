__author__ = 'hok1'

import sim_mr_canonical as sim
import csv

N = 10000000
total_energies = [10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 7500, 10000, 25000, 50000, 75000, 100000]
num_workers = 10
nrep = 2
filename = 'raw_canonical_sim.csv'

if __name__ == '__main__':
    simulator = sim.CanonicalEnsembleSimulationWrapper()
    outf = open(filename, 'wb')
    writer = csv.writer(outf)
    header = ['N', 'total_energy', 'mean_energy', 'std_energy', 'beta', 'alpha']
    writer.writerow(header)
    for total_energy in total_energies:
        print "Total Energy = ", total_energy
        for i in range(nrep):
            print "\tRepeat ", i
            occupancies = simulator.simulate_multithread(N, total_energy, numthreads=num_workers)
            coefs, mean_energy, std_energy = simulator.stat_analyze(occupancies)
            writer.writerow([N, total_energy, mean_energy, std_energy, -coefs[0], coefs[1]])
    outf.close()