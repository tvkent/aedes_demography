#!/bin/bash
#SBATCH -D /work/users/t/v/tvkent/aedes/Data/
#SBATCH --time=48:00:00
#SBATCH --mem=300G
#SBATCH -J admixture
#SBATCH --array=2-10
#SBATCH -o /work/users/t/v/tvkent/aedes/logs/admix-out-%A.%a.txt
#SBATCH -e /work/users/t/v/tvkent/aedes/logs/admix-error-%A.%a.txt

set -e
set -u

module load admixture/1.3.0


admixture \
	--cv \
	AaegL5_full_filtered_snps_110122_rep_map_masked_100kbgenes_centromere_filtered_100ld.bed \
	${SLURM_ARRAY_TASK_ID}


#-----------------------------------#
# Americas only
#-----------------------------------#


admixture \
	--cv \
	AaegL5_full_filtered_snps_110122_rep_map_masked_100kbgenes_centromere_filtered_100ld_americas.bed \
	${SLURM_ARRAY_TASK_ID}
