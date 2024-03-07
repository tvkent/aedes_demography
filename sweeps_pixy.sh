#!/bin/bash
#SBATCH -D /work/users/t/v/tvkent/aedes/Data/VCFs/
#SBATCH --time=1:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=4
#SBATCH -J pixy
#SBATCH --array=0-5
#SBATCH -o /work/users/t/v/tvkent/aedes/logs/pixy-out-%A.%a.txt
#SBATCH -e /work/users/t/v/tvkent/aedes/logs/pixy-error-%A.%a.txt

set -e
set -u

module load anaconda

conda activate pixy

module load bedtools


pop_id=${SLURM_ARRAY_TASK_ID}
#chr=$((pop_id%3))
var=$((SLURM_ARRAY_TASK_ID%2))

# n populations to calc pi on 
population=("Brazil" "Colombia" "USA" "Gabon" "Kenya" "Senegal") 
var=("0fold" "4fold")
# 6 vcfs to use
#vcfs=("AaegL5_1_invariant_variant_rep_map_depth_masked_0fold_fix.vcf.gz" "AaegL5_2_invariant_variant_rep_map_depth_masked_0fold_fix.vcf.gz" "AaegL5_3_invariant_variant_rep_map_depth_masked_0fold_fix.vcf.gz" "AaegL5_1_invariant_variant_rep_map_depth_masked_4fold_fix.vcf.gz" "AaegL5_2_invariant_variant_rep_map_depth_masked_4fold_fix.vcf.gz" "AaegL5_3_invariant_variant_rep_map_depth_masked_4fold_fix.vcf.gz")


python ../../Scripts/get_sweepfinder_bed.py ../Becca_sweepfinder/new_dataset/merged/AaegL5_1_${population[pop_id]}_merged_deduped_trimmed.121422.csv ../Becca_sweepfinder/new_dataset/merged/AaegL5_2_${population[pop_id]}_merged_deduped_trimmed.121422.csv ../Becca_sweepfinder/new_dataset/merged/AaegL5_3_${population[pop_id]}_merged_deduped_trimmed.121422.csv ../${population[pop_id]}_sweeps.01.bed

bedtools merge -i ../${population[pop_id]}_sweeps.01.bed > ../${population[pop_id]}_sweeps.01_merged.bed

pixy --stats pi \
--vcf AaegL5_full_invariant_variant_rep_map_depth_masked_fix.vcf.gz \
--populations ../populations.txt \
--bed_file ../${population[pop_id]}_sweeps.01_merged.bed \
--sites_file ../AaegL5_degenerate_0fold.txt \
--n_cores ${SLURM_CPUS_PER_TASK} \
--output_folder ../../Results/sweepfinder_diversity/ \
--output_prefix ${population[pop_id]}_0fold_sweeps.01

pixy --stats pi \
--vcf AaegL5_full_invariant_variant_rep_map_depth_masked_fix.vcf.gz \
--populations ../populations.txt \
--bed_file ../${population[pop_id]}_sweeps.01_merged.bed \
--sites_file ../AaegL5_degenerate_4fold.txt \
--n_cores ${SLURM_CPUS_PER_TASK} \
--output_folder ../../Results/sweepfinder_diversity/ \
--output_prefix ${population[pop_id]}_4fold_sweeps.01
