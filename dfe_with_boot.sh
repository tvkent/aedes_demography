#!/bin/bash
#SBATCH -D /work/users/t/v/tvkent/aedes/
#SBATCH --time=10:00:00
#SBATCH --mem=50G
#SBATCH -J dfe
#SBATCH --array=0-13
#SBATCH -o /work/users/t/v/tvkent/aedes/logs/dfe-out-%A.%a.txt
#SBATCH -e /work/users/t/v/tvkent/aedes/logs/dfe-error-%A.%a.txt
#--ntasks=24
#--array=0-13
set -e
set -u

module load anaconda

conda activate allel


populations=("FEMALE_1-F1_CGCATGAT-TCAGGCTT_S1,FEMALE_10-F10_GTGCCATA-ACTAGGAG_S2,FEMALE_11-F11_CGTTGCAA-CGCTCTAT_S3,FEMALE_12-F12_TGAAGACG-TGGCATGT_S4,FEMALE_14-F14_ACGTTCAG-GCACAACT_S6,FEMALE_15-F15_ATGCACGA-GACGATCT_S7,FEMALE_17-F17_CGGCTAAT-AGAACGAG_S9,FEMALE_18-F18_GAATCCGA-CACTAGCT_S10,FEMALE_2-F2_CTTAGGAC-GTAGGAGT_S12,FEMALE_20-F20_GTTACGCA-ATCGCCAT_S13,FEMALE_21-F21_ATGACGTC-GAAGGAAG_S14,FEMALE_22-F22_CAGTCCAA-GATTACCG_S15,FEMALE_3-F3_ATCCGGTA-TATCGGTC_S16,FEMALE_5-F5_CATACCAC-CGGTTGTT_S18,FEMALE_6-F6_AACGTCTG-GCGTCATT_S19,FEMALE_7-F7_CCTGATTG-AACTGAGC_S20,FEMALE_8-F8_CCTTGTAG-TTCCAAGG_S21,FEMALE_9-F9_ACGGAACA-GTTCTCGT_S22,MALE_1-M1_CGACGTTA-ATCCAGAG_S23,MALE_2-M2_TAACCGGT-GAGACGAT_S24,MALE_3-M3_ATCGATCG-TGCTTCCA_S25,MALE_4-M4_TCGCTGTT-ACGACTTG_S26,MALE_5-M5_CATGGCTA-GATGTGTG_S27,MALE_6-M6_AGCGTGTT-TTGCGAAG_S28" "JB_A2_18_S1,JB_A2_19_S2,JB_A2_25_S3,JB_A2_29_S4,JB_A2_34_S5,JB_A2_35_S6,JB_A2_36_S7,JB_A2_39_S8,JB_A2_46_S9,JB_A2_50_S10" "SRR11006835,SRR11006836,SRR11006837,SRR11006838,SRR11006839,SRR11006840,SRR11006841,SRR11006842,SRR11006843,SRR11006846,SRR11006847,SRR11006848,SRR11006849,SRR11006850,SRR11006851,SRR11006852,SRR11006853,SRR11006854" "SRR6768001,SRR6768002,SRR6768003,SRR6768004,SRR6768006,SRR6768007,SRR6768008,SRR6768009,SRR6768010,SRR6768011,SRR6768012,SRR6768013,SRR6768014,SRR6768015,SRR6768016,SRR6768017,SRR6768018,SRR6768019,SRR6768020,SRR6768021,SRR6768022,SRR6768023,SRR6768024,SRR6768025,SRR6768026,SRR6768027,SRR6768028" "SRR11006820,SRR11006821,SRR11006823,SRR11006824,SRR11006825,SRR11006826,SRR11006827,SRR11006828,SRR11006829,SRR11006830,SRR11006831,SRR11006832,SRR11006834" "SRR11006665,SRR11006666,SRR11006667,SRR11006668,SRR11006669,SRR11006670,SRR11006672,SRR11006673,SRR11006674,SRR11006675,SRR11006676,SRR11006677,SRR11006678,SRR11006679,SRR11006680,SRR11006681,SRR11006683,SRR11006684,SRR11006685" "SRR11006749,SRR11006751,SRR11006752,SRR11006753,SRR11006754,SRR11006755,SRR11006756,SRR11006757,SRR11006758,SRR11006759,SRR11006760,SRR11006762,SRR11006763,SRR11006764,SRR11006765,SRR11006766,SRR11006767,SRR11006768,SRR11006769,SRR11006770")

pop_names=('Rio_Claro' 'Cali' 'Brazil' 'USA' 'Gabon' 'Kenya' 'Senegal')

neut_site_set=('4fold' 'intergenic')
#neut_site_file=('masked_4fold' '10kb_intergenic_masked')
neut_site_file=('masked_4fold_fix' '100kbgenes_centromere_filtered_sites_mask')

site_ix=$((SLURM_ARRAY_TASK_ID % 2))
pop_ix=$((SLURM_ARRAY_TASK_ID / 2))
boot_ix=$((SLURM_ARRAY_TASK_ID / 2))


function run_dfe {

/nas/longleaf/home/tvkent/Software/dfe-alpha-release-2.15/est_dfe -c ${1}/boot${2}/neut_config_${2}.txt
/nas/longleaf/home/tvkent/Software/dfe-alpha-release-2.15/est_dfe -c ${1}/boot${2}/sel_config_${2}.txt

/nas/longleaf/home/tvkent/Software/dfe-alpha-release-2.15/prop_muts_in_s_ranges -c ${1}/boot${2}/sel/est_dfe.out -o ${1}/boot${2}/boot${2}.nes

}

export -f run_dfe


mkdir -p ./Results/DFE/${pop_names[pop_ix]}_${neut_site_file[site_ix]}

# 0: run dfe for real data
#whole genome
#AaegL5_full_filtered_snps_110122_rep_map_100kbgenes_centromere_filtered_sites_mask.vcf
python ./Scripts/sfs_for_dfe.py -s ./Data/VCFs/AaegL5_full_invariant_variant_rep_map_depth_masked_0fold_fix.vcf.gz -n ./Data/VCFs/AaegL5_full_invariant_variant_rep_map_depth_${neut_site_file[site_ix]}.vcf.gz -o ./Results/DFE/${pop_names[pop_ix]}_${neut_site_file[site_ix]} -p ${populations[pop_ix]}

python ./Scripts/sfs_for_dfe.py -s ./Data/VCFs/AaegL5_full_invariant_variant_rep_map_depth_masked_0fold_fix.vcf.gz -n ./Data/VCFs/AaegL5_full_filtered_snps_110122_rep_map_100kbgenes_centromere_filtered_sites_mask.vcf.gz -o ./Results/DFE/${pop_names[pop_ix]}_${neut_site_set[site_ix]} -p ${populations[pop_ix]}

/nas/longleaf/home/tvkent/Software/dfe-alpha-release-2.15/est_dfe -c ./Results/DFE/${pop_names[pop_ix]}_${neut_site_file[site_ix]}/neut_config.txt
/nas/longleaf/home/tvkent/Software/dfe-alpha-release-2.15/est_dfe -c ./Results/DFE/${pop_names[pop_ix]}_${neut_site_file[site_ix]}/sel_config.txt

/nas/longleaf/home/tvkent/Software/dfe-alpha-release-2.15/prop_muts_in_s_ranges -c ./Results/DFE/${pop_names[pop_ix]}_${neut_site_file[site_ix]}/sel/est_dfe.out -o ./Results/DFE/${pop_names[pop_ix]}_${neut_site_file[site_ix]}/dfe.nes


# bootstrapping
# 1: make sfs and configs
python ./Scripts/sfs_for_dfe.py -s ./Data/VCFs/AaegL5_full_invariant_variant_rep_map_depth_masked_0fold_fix.vcf.gz -n ./Data/VCFs/AaegL5_full_invariant_variant_rep_map_depth_${neut_site_file[site_ix]}_fix.vcf.gz -o ./Results/DFE/${pop_names[pop_ix]}_${neut_site_file[site_ix]} -p ${populations[pop_ix]} -b 1000

#2: run dfe
# run this part on muliple cores
parallel -j 24 run_dfe ::: "./Results/DFE/${pop_names[pop_ix]}_${neut_site_set[site_ix]}" ::: $(seq 0 999)


#3: process dfe results
echo -e "0-1\t1-10\t10-100\t100-inf" > ./Results/DFE/${pop_names[pop_ix]}_${neut_site_set[site_ix]}_bootstrapped1000.nes;for i in {0..999};do cat ./Results/DFE/${pop_names[pop_ix]}_${neut_site_set[site_ix]}/boot${i}/boot${i}.nes | cut -d" " -f3,6,9,12 | tr ' ' "\t" >> ./Results/DFE/${pop_names[pop_ix]}_${neut_site_set[site_ix]}_bootstrapped1000.nes;done

