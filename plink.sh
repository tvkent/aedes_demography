#!/bin/bash
#SBATCH -D /work/users/t/v/tvkent/aedes/Data/
#SBATCH -J plink
#SBATCH --time=10:00:00
#SBATCH --partition=dschridelab
#SBATCH --mem=20G
#SBATCH -o /work/users/t/v/tvkent/aedes/logs/out-plink-%A.%a.txt
#SBATCH -e /work/users/t/v/tvkent/aedes/logs/error-plink-%A.%a.txt

set -e
set -u

module load plink/1.90b3
#module load bedtools

#bedtools intersect -header -v -a ./Data/VCFs/AaegL5_full_filtered_snps_110122_rep_map_masked.vcf.gz -b ./Data/AaegL5_genes_100000_rep_map_centromere_filtered_sites_mask.bed.gz | \
#	sed "s/^AaegL5_//g" | \
#	./Data/VCFs/AaegL5_full_filtered_snps_110122_rep_map_masked_genes_100k_centromere_filtered.vcf.gz

#zcat AaegL5_${SLURM_ARRAY_TASK_ID}_filtered_snps_110122.vcf.gz | \
#sed "s/^AaegL5_//g" | \
#gzip > AaegL5_${SLURM_ARRAY_TASK_ID}_filtered_snps_110122_numericChr.vcf.gz
#zcat AaegL5_full_filtered_snps_110122_rep_map_10kgenic_masked.vcf.gz | \
#sed "s/^AaegL5_//g" | \
#gzip > AaegL5_full_filtered_snps_110122_rep_map_10kgenic_masked_numericChr.vcf.gz

#plink --make-bed \
#	--vcf VCFs/AaegL5_full_filtered_snps_110122_rep_map_100kbgenes_centromere_filtered_sites_mask_numericchr.vcf.gz \
#	--out AaegL5_full_filtered_snps_110122_rep_map_masked_100kbgenes_centromere_filtered \
#	--set-missing-var-ids @:#\$1,\$2 \
#	--double-id \
#	--allow-extra-chr \
#	--remove remove_medellin_plink.txt

#plink \
#--bfile AaegL5_full_filtered_snps_110122_rep_map_masked_100kbgenes_centromere_filtered \
#--out AaegL5_full_filtered_snps_110122_rep_map_masked_100kbgenes_centromere_filtered_100ld \
#--allow-extra-chr \
#--make-founders \
#--indep-pairwise 100 10 0.5

#make new bed with snps in low LD
plink \
--bfile AaegL5_full_filtered_snps_110122_rep_map_masked_100kbgenes_centromere_filtered \
--extract AaegL5_full_filtered_snps_110122_rep_map_masked_100kbgenes_centromere_filtered_100ld.prune.in \
--make-bed \
--out AaegL5_full_filtered_snps_110122_rep_map_masked_100kbgenes_centromere_filtered_100ld \
--allow-extra-chr

plink \
--bfile AaegL5_full_filtered_snps_110122_rep_map_masked_100kbgenes_centromere_filtered \
--extract AaegL5_full_filtered_snps_110122_rep_map_masked_100kbgenes_centromere_filtered_100ld.prune.in \
--make-bed \
--out AaegL5_full_filtered_snps_110122_rep_map_masked_100kbgenes_centromere_filtered_100ld_americas \
--allow-extra-chr \
--keep americas_samples.txt


#pca on thinned snps
plink \
--bfile AaegL5_full_filtered_snps_110122_rep_map_masked_100kbgenes_centromere_filtered \
--pca \
--allow-extra-chr \
--out ../Results/AaegL5_full_filtered_snps_110122_rep_map_masked_100kbgenes_centromere_filtered_100ld

#pca on americas only
plink \
--bfile AaegL5_full_filtered_snps_110122_rep_map_masked_100kbgenes_centromere_filtered_100ld_americas \
--pca \
--allow-extra-chr \
--out ../Results/AaegL5_full_filtered_snps_110122_rep_map_masked_100kbgenes_centromere_filtered_100ld_americas

#plink --make-bed \
#--vcf AaegL5_${SLURM_ARRAY_TASK_ID}_filtered_snps_110122_numericChr.vcf.gz \
#--out AaegL5_${SLURM_ARRAY_TASK_ID}_filtered_snps_110122 \
#--set-missing-var-ids @:#\$1,\$2 \
#--double-id \
#--allow-extra-chr
