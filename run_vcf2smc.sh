#!/bin/bash
#SBATCH -D /work/users/t/v/tvkent/aedes/
#SBATCH --time=50:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=2
#SBATCH -J vcf2smc
#SBATCH -o /work/users/t/v/tvkent/aedes/logs/vcf2smc-out-%A.%a.txt
#SBATCH -e /work/users/t/v/tvkent/aedes/logs/vcf2smc-error-%A.%a.txt

set -e
set -u

workdir=/work/users/t/v/tvkent/aedes

# source arguments from config
source $1

module load singularity

# convert distinguished individuals string into array
# and check which population it is 
sampleslist1="${pop1}:${samples1}"
sampleslist2="${pop2}:${samples2}"

IFS=',' read -r -a samples1array <<< "${samples1}"
IFS=',' read -r -a samples2array <<< "${samples2}"

IFS=',' read -r -a distinguished1 <<< "${dist1}"
IFS=',' read -r -a distinguished2 <<< "${dist2}"

distind1len=${#distinguished1[@]}


if [ $SLURM_ARRAY_TASK_ID -lt $distind1len ]; then
	distidx=${distinguished1[SLURM_ARRAY_TASK_ID]}
	distind=${samples1array[distidx]}
	focal_pop=$pop1

else
	index=$(( SLURM_ARRAY_TASK_ID - distind1len ))
	distidx=${distinguished2[index]}
	distind=${samples2array[distidx]}
	focal_pop=$pop2
fi


#lazily run all 3 chromosomes in same job
if [ $focal_pop == $pop1 ]; then
	singularity exec --bind /work/users/t/v/tvkent/aedes/ ./smc++ smc++ vcf2smc --cores ${SLURM_CPUS_PER_TASK} ${vcf} -d ${distind} ${distind} -m ./Data/AaegL5_genes_100000_rep_map_centromere_filtered_sites_mask.bed.gz ${outdir}/AaegL5_1_${pop1}-${pop2}_${distind}.smc.gz AaegL5_1 ${sampleslist1} ${sampleslist2}

	singularity exec --bind /work/users/t/v/tvkent/aedes/ ./smc++ smc++ vcf2smc --cores ${SLURM_CPUS_PER_TASK} ${vcf} -d ${distind} ${distind} -m ./Data/AaegL5_genes_100000_rep_map_centromere_filtered_sites_mask.bed.gz ${outdir}/AaegL5_2_${pop1}-${pop2}_${distind}.smc.gz AaegL5_2 ${sampleslist1} ${sampleslist2}

	singularity exec --bind /work/users/t/v/tvkent/aedes/ ./smc++ smc++ vcf2smc --cores ${SLURM_CPUS_PER_TASK} ${vcf} -d ${distind} ${distind} -m ./Data/AaegL5_genes_100000_rep_map_centromere_filtered_sites_mask.bed.gz ${outdir}/AaegL5_3_${pop1}-${pop2}_${distind}.smc.gz AaegL5_3 ${sampleslist1} ${sampleslist2}

else
	singularity exec --bind /work/users/t/v/tvkent/aedes/ ./smc++ smc++ vcf2smc --cores ${SLURM_CPUS_PER_TASK} ${vcf} -d ${distind} ${distind} -m ./Data/AaegL5_genes_100000_rep_map_centromere_filtered_sites_mask.bed.gz ${outdir}/AaegL5_1_${pop2}-${pop1}_${distind}.smc.gz AaegL5_1 ${sampleslist2} ${sampleslist1}

	singularity exec --bind /work/users/t/v/tvkent/aedes/ ./smc++ smc++ vcf2smc --cores ${SLURM_CPUS_PER_TASK} ${vcf} -d ${distind} ${distind} -m ./Data/AaegL5_genes_100000_rep_map_centromere_filtered_sites_mask.bed.gz ${outdir}/AaegL5_2_${pop2}-${pop1}_${distind}.smc.gz AaegL5_2 ${sampleslist2} ${sampleslist1}

	singularity exec --bind /work/users/t/v/tvkent/aedes/ ./smc++ smc++ vcf2smc --cores ${SLURM_CPUS_PER_TASK} ${vcf} -d ${distind} ${distind} -m ./Data/AaegL5_genes_100000_rep_map_centromere_filtered_sites_mask.bed.gz ${outdir}/AaegL5_3_${pop2}-${pop1}_${distind}.smc.gz AaegL5_3 ${sampleslist2} ${sampleslist1}

fi
