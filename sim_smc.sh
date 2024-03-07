#!/bin/bash
#SBATCH -D /work/users/t/v/tvkent/aedes/
#SBATCH --time=5:00:00
#SBATCH --mem=5G
#SBATCH -J sfs_sim
#SBATCH --array=6
#SBATCH -o /work/users/t/v/tvkent/aedes/logs/sfs-sim-out-%A.%a.txt
#SBATCH -e /work/users/t/v/tvkent/aedes/logs/sfs-sim-error-%A.%a.txt

set -e
set -u

module load anaconda

conda activate ts

pop_names=("rio_claro" "cali" "kenya" "senegal" "gabon" "brazil" "us")
smc_names=("Rio_Claro" "Cali" "Kenya" "Senegal" "Gabon" "Brazil" "USA")

python3 Scripts/sim_smc.py \
	-s Results/${pop_names[SLURM_ARRAY_TASK_ID]}_rep_map_100kbgenes_centromere_filtered_sites_mask_folded.sfs \
	-d Results/smc++/${smc_names[SLURM_ARRAY_TASK_ID]}_100k_sites/${smc_names[SLURM_ARRAY_TASK_ID]}_100ksites_smc.csv \
	-o Results/smc++/${smc_names[SLURM_ARRAY_TASK_ID]}_100k_sites_rp5_eval.txt

python3 Scripts/sim_smc.py \
        -s Results/${pop_names[SLURM_ARRAY_TASK_ID]}_rep_map_100kbgenes_centromere_filtered_sites_mask_folded.sfs \
        -d Results/smc++/${smc_names[SLURM_ARRAY_TASK_ID]}_100k_sites_rp4/${smc_names[SLURM_ARRAY_TASK_ID]}_100ksites_rp4_smc.csv \
        -o Results/smc++/${smc_names[SLURM_ARRAY_TASK_ID]}_100k_sites_rp4_eval.txt

python3 Scripts/sim_smc.py \
        -s Results/${pop_names[SLURM_ARRAY_TASK_ID]}_rep_map_100kbgenes_centromere_filtered_sites_mask_folded.sfs \
        -d Results/smc++/${smc_names[SLURM_ARRAY_TASK_ID]}_100k_sites_rp6/${smc_names[SLURM_ARRAY_TASK_ID]}_100ksites_rp6_smc.csv \
        -o Results/smc++/${smc_names[SLURM_ARRAY_TASK_ID]}_100k_sites_rp6_eval.txt

