#!/bin/bash
#SBATCH -D /work/users/t/v/tvkent/aedes/Data/VCFs/
#SBATCH --time=20:00:00
#SBATCH --mem=3G
#SBATCH --cpus-per-task=4
#SBATCH -J pixy
#SBATCH --array=0-2
#SBATCH -o /work/users/t/v/tvkent/aedes/logs/pixy-out-%A.%a.txt
#SBATCH -e /work/users/t/v/tvkent/aedes/logs/pixy-error-%A.%a.txt

set -e
set -u

module load anaconda

conda activate pixy

#     ${SLURM_ARRAY_TASK_ID}

#vcfs=("AaegL5_1_invariant_variant_rep_map_depth_masked_0fold_fix.vcf.gz" "AaegL5_2_invariant_variant_rep_map_depth_masked_0fold_fix.vcf.gz" "AaegL5_3_invariant_variant_rep_map_depth_masked_0fold_fix.vcf.gz" "AaegL5_1_invariant_variant_rep_map_depth_masked_4fold_fix.vcf.gz" "AaegL5_2_invariant_variant_rep_map_depth_masked_4fold_fix.vcf.gz" "AaegL5_3_invariant_variant_rep_map_depth_masked_4fold_fix.vcf.gz" "AaegL5_1_invariant_variant_rep_map_depth_10kb_intergenic_masked_fix.vcf.gz" "AaegL5_2_invariant_variant_rep_map_depth_10kb_intergenic_masked_fix.vcf.gz" "AaegL5_3_invariant_variant_rep_map_depth_10kb_intergenic_masked_fix.vcf.gz" "AaegL5_1_invariant_variant_rep_map_depth_masked_fix.vcf.gz" "AaegL5_2_invariant_variant_rep_map_depth_masked_fix.vcf.gz" "AaegL5_3_invariant_variant_rep_map_depth_masked_fix.vcf.gz")
#vcfs=(0fold 4fold intergenic all)
vcfs=("AaegL5_1_invariant_variant_rep_map_depth_masked_fix.vcf.gz" "AaegL5_2_invariant_variant_rep_map_depth_masked_fix.vcf.gz" "AaegL5_3_invariant_variant_rep_map_depth_masked_fix.vcf.gz")


pixy --stats pi dxy fst \
--vcf ${vcfs[SLURM_ARRAY_TASK_ID]} \
--populations ../subpopulations.txt \
--n_cores ${SLURM_CPUS_PER_TASK} \
--output_folder ../../Results/ \
--output_prefix ${vcfs[SLURM_ARRAY_TASK_ID]}_wholegenome
