#!/bin/bash
#SBATCH -D /work/users/t/v/tvkent/aedes/Data/
#SBATCH --time=4:00:00
#SBATCH --mem=50G
#SBATCH -J mask
#SBATCH -o /work/users/t/v/tvkent/aedes/logs/mask-out-%A.%a.txt
#SBATCH -e /work/users/t/v/tvkent/aedes/logs/mask-error-%A.%a.txt

set -e
set -u

module load bedtools

gff=/proj/dschridelab/tvkent/aedes/Data/AaegL5.0_genomic_NCBI_renamed_genes.gff.gz
flank_size=100000

#get mask of genes and flank
cat <(bedtools flank -b ${flank_size} -i <(zgrep "^AaegL5" ${gff} | sort -k1,1 -k4,4n) -g ./VCFs/AaegL5_genome.txt | awk 'BEGIN{OFS=FS="\t"} {print $1,$4-1,$5}') <(zgrep "^AaegL5" ${gff} | sort -k1,1 -k4,4n | awk 'BEGIN{OFS=FS="\t"} {print $1,$4-1,$5}') | \
	sort -k1,1 -k2,2n | \
       	bedtools merge -i stdin > \
       	AaegL5_genes_${flank_size}_mask.bed

#merge with becca's mask and centromere locations
cat AaegL5_genes_${flank_size}_mask.bed /overflow/dschridelab/users/rrlove/aedes/refs/aegy/unified_mask/merged_rep_map_masks.110822.bed AaegL5_centromeres.bed |\
       	sort -k1,1 -k2,2n | \
	bedtools merge -i stdin > \
       	AaegL5_genes_${flank_size}_rep_map_centromere_mask.bed

#merge with complement of filtered snps and invariants
#should keep only snps and quality ROH, everything else as missing
cat <(bedtools complement -i ./VCFs/AaegL5_full_invariant_variant_rep_map_depth_masked_fix.vcf.gz -g ./VCFs/AaegL5_genome.txt) AaegL5_genes_${flank_size}_rep_map_centromere_mask.bed | \
	sort -k1,1 -k2,2n | \
	bedtools merge -i stdin > \
	AaegL5_genes_${flank_size}_rep_map_centromere_filtered_sites_mask.bed
