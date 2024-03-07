#!/bin/bash
#SBATCH -D /proj/dschridelab/tvkent/aedes/
#SBATCH --time=10:00:00
#SBATCH --mem=150G
#SBATCH -J concat
#SBATCH -o /proj/dschridelab/tvkent/aedes/out-%A.%a.txt
#SBATCH -e /proj/dschridelab/tvkent/aedes/error-%A.%a.txt

set -e
set -u

module use $HOME/modulefiles
module load bcftools

bcftools concat -a -Oz -o AaegL5_full_filtered_snps_110122.vcf.gz AaegL5_1_filtered_snps_110122.vcf.gz AaegL5_2_filtered_snps_110122.vcf.gz AaegL5_3_filtered_snps_110122.vcf.gz
