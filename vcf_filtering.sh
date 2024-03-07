#!/bin/bash
#SBATCH -D /work/users/t/v/tvkent/aedes/
#SBATCH --time=40:00:00
#SBATCH --mem=50G
#SBATCH -J degen
#SBATCH -o /work/users/t/v/tvkent/aedes/logs/fix_lg3-out-%A.%a.txt
#SBATCH -e /work/users/t/v/tvkent/aedes/logs/fix_lg3-error-%A.%a.txt
set -e
set -u

module use $HOME/modulefiles
module load bcftools
module load bedtools/2.30
module load vcftools/0.1.15

projdir=/proj/dschridelab/tvkent/aedes

#filter out repetitive and poorly-mapped regions
#bedtools intersect -v -header -wa -a ${projdir}/AaegL5_${SLURM_ARRAY_TASK_ID}_filtered_snps_110122.vcf.gz -b /overflow/dschridelab/users/rrlove/aedes/refs/aegy/unified_mask/merged_rep_map_masks.110822.bed | bgzip -c > ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_filtered_snps_110122_rep_map_masked.vcf.gz

#tabix -p vcf ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_filtered_snps_110122_rep_map_masked.vcf.gz

#filter out 10kb around genes
#bedtools intersect -header -wa -a ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_filtered_snps_110122_rep_map_masked.vcf.gz -b ${projdir}/Data/AaegL5_10kb_intergenic.bed | bgzip -c > ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_filtered_snps_110122_rep_map_10kgenic_masked.vcf.gz

#tabix -p vcf ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_filtered_snps_110122_rep_map_10kgenic_masked.vcf.gz

#merge and index
#run this on one core
#bcftools concat --naive -O z -o ./Data/VCFs/AaegL5_full_filtered_snps_110122_rep_map_masked.vcf.gz ./Data/VCFs/AaegL5_1_filtered_snps_110122_rep_map_masked.vcf.gz ./Data/VCFs/AaegL5_2_filtered_snps_110122_rep_map_masked.vcf.gz ./Data/VCFs/AaegL5_3_filtered_snps_110122_rep_map_masked.vcf.gz

#tabix -p vcf ./Data/VCFs/AaegL5_full_filtered_snps_110122_rep_map_masked.vcf.gz

#bcftools concat --naive -O z -o ./Data/VCFs/AaegL5_full_filtered_snps_110122_rep_map_10kgenic_masked.vcf.gz ./Data/VCFs/AaegL5_1_filtered_snps_110122_rep_map_10kgenic_masked.vcf.gz ./Data/VCFs/AaegL5_2_filtered_snps_110122_rep_map_10kgenic_masked.vcf.gz ./Data/VCFs/AaegL5_3_filtered_snps_110122_rep_map_10kgenic_masked.vcf.gz

#tabix -p vcf ./Data/VCFs/AaegL5_full_filtered_snps_110122_rep_map_10kgenic_masked.vcf.gz

#merge variant and invariants
#back to multicore
#tabix -p vcf ${projdir}/AaegL5_full_filtered_snps_110122.vcf.gz

#split invariants by chromosome and remove sites with more than 25% samples missing
#bcftools view -O b -r AaegL5_${SLURM_ARRAY_TASK_ID} -e 'F_MISSING>0.25' whole_sample_all_confident_AaegL5_full_invariant_frombcf_nolowqual.vcf.gz > ./Data/VCFs/whole_sample_all_confident_AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_frombcf_nolowqual.bcf

#bcftools index ./Data/VCFs/whole_sample_all_confident_AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_frombcf_nolowqual.bcf

#bcftools concat -a -O z -o ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant.vcf.gz ${projdir}/AaegL5_full_filtered_snps_110122.vcf.gz ./Data/VCFs/whole_sample_all_confident_AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_frombcf_nolowqual.bcf

#tabix -p vcf ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant.vcf.gz
#bcftools view -r AaegL5_${SLURM_ARRAY_TASK_ID} -O z -o ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_fix.vcf.gz ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant.vcf.gz

#tabix -p vcf ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_fix.vcf.gz

#bcftools concat --naive -O z -o ./Data/VCFs/AaegL5_full_invariant_variant_fix.vcf.gz ./Data/VCFs/AaegL5_1_invariant_variant_fix.vcf.gz ./Data/VCFs/AaegL5_2_invariant_variant_fix.vcf.gz ./Data/VCFs/AaegL5_3_invariant_variant_fix.vcf.gz

#tabix -p vcf ./Data/VCFs/AaegL5_full_invariant_variant_fix.vcf.gz
#remove repetitive and poorly-mapped regions
#bedtools intersect -v -header -wa -a ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant.vcf.gz -b /overflow/dschridelab/users/rrlove/aedes/refs/aegy/unified_mask/merged_rep_map_masks.110822.bed | bgzip -c > ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_masked.vcf.gz

#tabix -p vcf ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_masked.vcf.gz
#bcftools view -r AaegL5_${SLURM_ARRAY_TASK_ID} -O z -o ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_masked_fix.vcf.gz ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_masked.vcf.gz

#tabix -p vcf ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_masked_fix.vcf.gz

#bcftools concat --naive -O z -o ./Data/VCFs/AaegL5_full_invariant_variant_rep_map_masked_fix.vcf.gz ./Data/VCFs/AaegL5_1_invariant_variant_rep_map_masked_fix.vcf.gz ./Data/VCFs/AaegL5_2_invariant_variant_rep_map_masked_fix.vcf.gz ./Data/VCFs/AaegL5_3_invariant_variant_rep_map_masked_fix.vcf.gz

#tabix -p vcf ./Data/VCFs/AaegL5_full_invariant_variant_rep_map_masked_fix.vcf.gz
#remove low depth sites
#bedtools intersect -v -header -wa -a ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_masked.vcf.gz -b ./Data/whole_sample_all_confident_AaegL5_1_invariant_frombcf_nolowqual_missingfilt_depthmask_4-30.bed | bgzip -c > ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_masked.vcf.gz

#tabix -p vcf ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_masked.vcf.gz
#bcftools view -r AaegL5_${SLURM_ARRAY_TASK_ID} -O z -o ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_masked_fix.vcf.gz ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_masked.vcf.gz

#tabix -p vcf ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_masked_fix.vcf.gz

#bcftools concat --naive -O z -o ./Data/VCFs/AaegL5_full_invariant_variant_rep_map_depth_masked_fix.vcf.gz ./Data/VCFs/AaegL5_1_invariant_variant_rep_map_depth_masked_fix.vcf.gz ./Data/VCFs/AaegL5_2_invariant_variant_rep_map_depth_masked_fix.vcf.gz ./Data/VCFs/AaegL5_3_invariant_variant_rep_map_depth_masked_fix.vcf.gz

#tabix -p vcf ./Data/VCFs/AaegL5_full_invariant_variant_rep_map_depth_masked_fix.vcf.gz
#remove 10kb around genes
#bedtools intersect -v -header -wa -a ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_masked.vcf.gz -b ${projdir}/Data/AaegL5_10kb_intergenic.bed | bgzip -c > ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_10kb_intergenic_masked.vcf.gz

#tabix -p vcf ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_10kb_intergenic_masked.vcf.gz
#bcftools view -r AaegL5_${SLURM_ARRAY_TASK_ID} -O z -o ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_10kb_intergenic_masked_fix.vcf.gz ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_10kb_intergenic_masked.vcf.gz

#tabix -p vcf ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_10kb_intergenic_masked_fix.vcf.gz

#bcftools concat --naive -O z -o ./Data/VCFs/AaegL5_full_invariant_variant_rep_map_depth_10kb_intergenic_masked_fix.vcf.gz ./Data/VCFs/AaegL5_1_invariant_variant_rep_map_depth_10kb_intergenic_masked_fix.vcf.gz ./Data/VCFs/AaegL5_2_invariant_variant_rep_map_depth_10kb_intergenic_masked_fix.vcf.gz ./Data/VCFs/AaegL5_3_invariant_variant_rep_map_depth_10kb_intergenic_masked_fix.vcf.gz

#tabix -p vcf ./Data/VCFs/AaegL5_full_invariant_variant_rep_map_depth_10kb_intergenic_masked_fix.vcf.gz

# 4fold and 0fold sites
#bedtools intersect -header -wa -a ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_masked.vcf.gz -b ${projdir}/Data/AaegL5_degenerate_full_0fold.bed | bgzip -c > ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_masked_0fold.vcf.gz

#tabix -p vcf ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_masked_0fold.vcf.gz

#bcftools view -r AaegL5_${SLURM_ARRAY_TASK_ID} -O z -o ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_masked_0fold_fix.vcf.gz ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_masked_0fold.vcf.gz

#tabix -p vcf ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_masked_0fold_fix.vcf.gz

bcftools concat --naive -O z -o ./Data/VCFs/AaegL5_full_invariant_variant_rep_map_depth_masked_0fold_fix.vcf.gz ./Data/VCFs/AaegL5_1_invariant_variant_rep_map_depth_masked_0fold_fix.vcf.gz ./Data/VCFs/AaegL5_2_invariant_variant_rep_map_depth_masked_0fold_fix.vcf.gz ./Data/VCFs/AaegL5_3_invariant_variant_rep_map_depth_masked_0fold_fix.vcf.gz

tabix -p vcf ./Data/VCFs/AaegL5_full_invariant_variant_rep_map_depth_masked_0fold_fix.vcf.gz
#bedtools intersect -header -wa -a ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_masked.vcf.gz -b ${projdir}/Data/AaegL5_degenerate_full_4fold.bed | bgzip -c > ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_masked_4fold.vcf.gz

#tabix -p vcf ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_masked_4fold.vcf.gz

#bcftools view -r AaegL5_${SLURM_ARRAY_TASK_ID} -O z -o ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_masked_4fold_fix.vcf.gz ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_masked_4fold.vcf.gz

#tabix -p vcf ./Data/VCFs/AaegL5_${SLURM_ARRAY_TASK_ID}_invariant_variant_rep_map_depth_masked_4fold_fix.vcf.gz

bcftools concat --naive -O z -o ./Data/VCFs/AaegL5_full_invariant_variant_rep_map_depth_masked_4fold_fix.vcf.gz ./Data/VCFs/AaegL5_1_invariant_variant_rep_map_depth_masked_4fold_fix.vcf.gz ./Data/VCFs/AaegL5_2_invariant_variant_rep_map_depth_masked_4fold_fix.vcf.gz ./Data/VCFs/AaegL5_3_invariant_variant_rep_map_depth_masked_4fold_fix.vcf.gz

tabix -p vcf ./Data/VCFs/AaegL5_full_invariant_variant_rep_map_depth_masked_4fold_fix.vcf.gz
