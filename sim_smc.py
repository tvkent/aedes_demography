import tskit
import msprime
import pandas as pd
import numpy as np
import math
import argparse

def main_arguments():
    parser = argparse.ArgumentParser(description="Simulate smc++ demography and compare SFS with LL")
    parser.add_argument('-s', '--obssfs', help='observed sfs',required=True)
    parser.add_argument('-d', '--smc', help='smc++ csv',required=True)
    parser.add_argument('-o', '--output', help='output path',required=True)
    args = parser.parse_args()
    return(args)

args = main_arguments()

# read observed folded sfs
obs_sfs = pd.read_csv(args.obssfs, sep=" ",header=None).values[0]
obs_sfs=np.delete(obs_sfs,0)
nsamples=obs_sfs.shape[0]

# read in smc csv for demog
size_history = pd.read_csv(args.smc)

sizes = size_history['y']
times = size_history['x'][1:]

# define demography as smc results
demography=msprime.Demography()
demography.add_population(name="A", initial_size=sizes[0])

for pop_size, time in zip(sizes[1:], times):
    demography.add_population_parameters_change(time=time*15,initial_size=pop_size,population='A')

# simulate segment under smc demography
tree_sequence = msprime.sim_ancestry(sequence_length=100000,recombination_rate=6e-7,demography=demography,samples=nsamples)
mts = msprime.sim_mutations(tree_sequence, rate=4.85e-9)

af = mts.allele_frequency_spectrum()

# get rid of folded zeros and rescale as freq
sim_sfs=af[af != 0]
sim_freq=sim_sfs/np.sum(sim_sfs)

# get residuals from frequency sfs
obs_freq=obs_sfs/np.sum(obs_sfs)

resid=np.subtract(obs_freq,sim_freq)
# get multinomial log likelihood of simulated sfs
llSum_sim=0
for i in range(0,obs_sfs.shape[0]):
    llSum_sim += obs_sfs[i]*math.log(sim_freq[i])

llSum_self=0
for i in range(0,obs_sfs.shape[0]):
    llSum_self += obs_sfs[i]*math.log(obs_freq[i])

with open(args.output, "w") as f:
    f.write(str(llSum_sim) + "\n" + str(llSum_self) + "\n")

with open(args.output, "a") as f:
    np.savetxt(f, sim_sfs, newline="\t")
    f.write("\n")

with open(args.output, "a") as f:
    np.savetxt(f, obs_sfs, newline="\t")
    f.write("\n")
    
with open(args.output, "a") as f:
    np.savetxt(f, resid, newline="\t")


