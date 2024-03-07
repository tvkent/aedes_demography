#!/bin/bash
#SBATCH -D /work/users/t/v/tvkent/aedes/
#SBATCH --time=15:00:00
#SBATCH --mem=100G
#SBATCH -J stairway
#SBATCH --array=0-8
#SBATCH -o /work/users/t/v/tvkent/aedes/logs/stairway-out-%A.%a.txt
#SBATCH -e /work/users/t/v/tvkent/aedes/logs/stairway-error-%A.%a.txt

set -e
set -u

workdir=/work/users/t/v/tvkent/aedes

rio_claro="FEMALE_1-F1_CGCATGAT-TCAGGCTT_S1,FEMALE_10-F10_GTGCCATA-ACTAGGAG_S2,FEMALE_11-F11_CGTTGCAA-CGCTCTAT_S3,FEMALE_12-F12_TGAAGACG-TGGCATGT_S4,FEMALE_14-F14_ACGTTCAG-GCACAACT_S6,FEMALE_15-F15_ATGCACGA-GACGATCT_S7,FEMALE_17-F17_CGGCTAAT-AGAACGAG_S9,FEMALE_18-F18_GAATCCGA-CACTAGCT_S10,FEMALE_2-F2_CTTAGGAC-GTAGGAGT_S12,FEMALE_20-F20_GTTACGCA-ATCGCCAT_S13,FEMALE_21-F21_ATGACGTC-GAAGGAAG_S14,FEMALE_22-F22_CAGTCCAA-GATTACCG_S15,FEMALE_3-F3_ATCCGGTA-TATCGGTC_S16,FEMALE_5-F5_CATACCAC-CGGTTGTT_S18,FEMALE_6-F6_AACGTCTG-GCGTCATT_S19,FEMALE_7-F7_CCTGATTG-AACTGAGC_S20,FEMALE_8-F8_CCTTGTAG-TTCCAAGG_S21,FEMALE_9-F9_ACGGAACA-GTTCTCGT_S22,MALE_1-M1_CGACGTTA-ATCCAGAG_S23,MALE_2-M2_TAACCGGT-GAGACGAT_S24,MALE_3-M3_ATCGATCG-TGCTTCCA_S25,MALE_4-M4_TCGCTGTT-ACGACTTG_S26,MALE_5-M5_CATGGCTA-GATGTGTG_S27,MALE_6-M6_AGCGTGTT-TTGCGAAG_S28"

cali="JB_A2_18_S1,JB_A2_19_S2,JB_A2_25_S3,JB_A2_29_S4,JB_A2_34_S5,JB_A2_35_S6,JB_A2_36_S7,JB_A2_39_S8,JB_A2_46_S9,JB_A2_50_S10"

kenya="SRR11006665,SRR11006666,SRR11006667,SRR11006668,SRR11006669,SRR11006670,SRR11006672,SRR11006673,SRR11006674,SRR11006675,SRR11006676,SRR11006677,SRR11006678,SRR11006679,SRR11006680,SRR11006681,SRR11006683,SRR11006684,SRR11006685"

senegal="SRR11006749,SRR11006751,SRR11006752,SRR11006753,SRR11006754,SRR11006755,SRR11006756,SRR11006757,SRR11006758,SRR11006759,SRR11006760,SRR11006762,SRR11006763,SRR11006764,SRR11006765,SRR11006766,SRR11006767,SRR11006768,SRR11006769,SRR11006770"

gabon="SRR11006820,SRR11006821,SRR11006823,SRR11006824,SRR11006825,SRR11006826,SRR11006827,SRR11006828,SRR11006829,SRR11006830,SRR11006831,SRR11006832,SRR11006834"

brazil="SRR11006835,SRR11006836,SRR11006837,SRR11006838,SRR11006839,SRR11006840,SRR11006841,SRR11006842,SRR11006843,SRR11006846,SRR11006847,SRR11006848,SRR11006849,SRR11006850,SRR11006851,SRR11006852,SRR11006853,SRR11006854"

clovis="SRR6768019,SRR6768020,SRR6768021,SRR6768022,SRR6768023,SRR6768025,SRR6768026,SRR6768028"

rest_of_us="SRR6768001,SRR6768002,SRR6768003,SRR6768004,SRR6768006,SRR6768007,SRR6768008,SRR6768009,SRR6768010,SRR6768011,SRR6768012,SRR6768013,SRR6768014,SRR6768015,SRR6768016,SRR6768017,SRR6768018,SRR6768024,SRR6768027"

us="SRR6768019,SRR6768020,SRR6768021,SRR6768022,SRR6768023,SRR6768025,SRR6768026,SRR6768028,SRR6768001,SRR6768002,SRR6768003,SRR6768004,SRR6768006,SRR6768007,SRR6768008,SRR6768009,SRR6768010,SRR6768011,SRR6768012,SRR6768013,SRR6768014,SRR6768015,SRR6768016,SRR6768017,SRR6768018,SRR6768024,SRR6768027"

populations=($rio_claro $cali $kenya $senegal $gabon $brazil $clovis $rest_of_us $us)
pop_names=("rio_claro" "cali" "kenya" "senegal" "gabon" "brazil" "clovis" "rest_of_us" "us")
nseq_pops=("24" "10" "19" "20" "13" "18" "8" "19" "27")

cd $HOME/Software/stairway_plot/v2.1.1/

java -cp stairway_plot_es Stairbuilder ${workdir}/Data/${pop_names[SLURM_ARRAY_TASK_ID]}_rep_map_10k_masked_folded.blueprint

bash ${workdir}/Data/${pop_names[SLURM_ARRAY_TASK_ID]}_rep_map_10k_masked_folded.blueprint.sh

