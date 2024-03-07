#!/bin/bash
#SBATCH -D /work/users/t/v/tvkent/aedes/
#SBATCH --time=30:00:00
#SBATCH --mem=4M
#SBATCH --array=0-1
#SBATCH --partition=dschridelab
#SBATCH -J smc++split
#SBATCH -o /work/users/t/v/tvkent/aedes/logs/smc++auto-out-%A.%a.txt
#SBATCH -e /work/users/t/v/tvkent/aedes/logs/smc++auto-error-%A.%a.txt
#0-34
set -e
set -u

workdir=/work/users/t/v/tvkent/aedes
subpopulations=true

#pop_pairs=("NorCal-Rio_Claro" "NorCal-Cali" "NorCal-Brazil" "NorCal-Senegal" "NorCal-Kenya" "NorCal-Gabon" "SoCal-Rio_Claro" "SoCal-Cali" "SoCal-Brazil" "Socal-Senegal" "SoCal-Kenya" "Socal-Gabon" "Florida-Rio_Claro" "Florida-Cali" "Florida-Brazil" "Florida-Senegal" "Florida-Kenya" "Florida-Gabon" "Clovis-Rio_Claro" "Clovis-Cali" "Clovis-Brazil" "Clovis-Senegal" "Clovis-Kenya" "Clovis-Gabon" "Exeter-Rio_Claro" "Exeter-Cali" "Exeter-Brazil" "Exeter-Senegal" "Exeter-Kenya" "Exeter-Gabon")

#pop_pairs=("SoCal-Rio_Claro" "NorCal-Rio_Claro" "Florida-Rio_Claro" "Exeter-Rio_Claro" "Clovis-Rio_Claro" "SoCal-Senegal" "SoCal-Gabon")
pop_pairs=("SoCal-Senegal" "SoCal-Gabon")
pop_pair=${pop_pairs[SLURM_ARRAY_TASK_ID]}

pop1=${pop_pair%-*}
pop2=${pop_pair#*-}

# search populations file for each pop to get samples
if [ "$subpopulations" = true ]; then

	# get list of samples as a comma-sep string at subpopulation level
	samples1=$( awk -v var=$pop1 'BEGIN{ORS=","} $3==var{print $1}' ./Data/populations_with_subpops.txt | rev | cut -c 2- | rev )
	samples2=$( awk -v var=$pop2 'BEGIN{ORS=","} $3==var{print $1}' ./Data/populations_with_subpops.txt | rev | cut -c 2- | rev )

else
	# get list of samples as comma-sep string at population level (USA as country)
	samples1=$( awk -v var=$pop1 'BEGIN{ORS=","} $2==var{print $1}' ./Data/populations_with_subpops.txt | rev | cut -c 2- | rev )
	samples2=$( awk -v var=$pop2 'BEGIN{ORS=","} $2==var{print $1}' ./Data/populations_with_subpops.txt | rev | cut -c 2- | rev )
fi

# string format for smc++ list of samples
sampleslist1="${pop1}:${samples1}"
sampleslist2="${pop2}:${samples2}"

# array of sample names
IFS=',' read -r -a samples1array <<< "${samples1}"
IFS=',' read -r -a samples2array <<< "${samples2}"

# length of array in 1 and 0-based 
array1len=${#samples1array[@]}
rand1=$(( array1len - 1 ))
array2len=${#samples2array[@]}
rand2=$(( array2len - 1 ))
echo $array1len
echo $array2len
# get array of distinguished invididual id's 
# 5 distinguished if more than 5 individuals
# 2 distinguished if 3 or 4
# 1 distinguished if 1 or 2
if [ $array1len -gt 5 ];then
	distidx1=$( shuf -i0-$rand1 -n5 | tr '\n' ',' | rev | cut -c 2- | rev )
	IFS=',' read -r -a distidx1array <<< "${distidx1}"

else
	distlen1=$(( rand1 - 1 ))
	if [ $distlen1 -ge 2 ]; then
		distidx1=$( shuf -i0-$rand1 -n2 | tr '\n' ',' | rev | cut -c 2- | rev )
		IFS=',' read -r -a distidx1array <<< "${distidx1}"

	else
		distidx1=$( shuf -i0-$rand1 -n1 | tr '\n' ',' | rev | cut -c 2- | rev )
		IFS=',' read -r -a distidx1array <<< "${distidx1}"
	fi
fi

if [ $array2len -gt 5 ];then
	distidx2=$( shuf -i0-$rand2 -n5 | tr '\n' ',' | rev | cut -c 2- | rev )
	IFS=',' read -r -a distidx2array <<< "${distidx2}"

else
	distlen2=$(( rand2 - 1 ))
	if [ $distlen2 -ge 2 ]; then
		distidx2=$( shuf -i0-$rand2 -n2 | tr '\n' ',' | rev | cut -c 2- | rev )
		IFS=',' read -r -a distidx2array <<< "${distidx2}"

	else
		distidx2=$( shuf -i0-$rand2 -n1 | tr '\n' ',' | rev | cut -c 2- | rev )
		IFS=',' read -r -a distidx2array <<< "${distidx2}"
	fi
fi

echo $distidx1
echo $distidx2

# run vcf2smc in parallel for all distinguished individuals
distind1len=${#distidx1array[@]}
distind2len=${#distidx2array[@]}
num=$(( distind1len + distind2len ))
num_array=$(( num - 1 ))
echo ${distidx1array[@]}
echo ${distidx2array[@]}
# make config file
mkdir -p Data/smc++/${pop_pair}/

echo "vcf=Data/VCFs/AaegL5_full_filtered_snps_110122_rep_map_masked.vcf.gz" > Data/smc++/${pop_pair}/${pop_pair}_vcf2smc.config
echo "mask=Data/AaegL5_genes_100000_rep_map_centromere_filtered_sites_mask.bed.gz" >> Data/smc++/${pop_pair}/${pop_pair}_vcf2smc.config
echo "samples1=${samples1}" >> Data/smc++/${pop_pair}/${pop_pair}_vcf2smc.config
echo "samples2=${samples2}" >> Data/smc++/${pop_pair}/${pop_pair}_vcf2smc.config
echo "pop=${pop_pair}" >> Data/smc++/${pop_pair}/${pop_pair}_vcf2smc.config
echo "pop1=${pop1}" >> Data/smc++/${pop_pair}/${pop_pair}_vcf2smc.config
echo "pop2=${pop2}" >> Data/smc++/${pop_pair}/${pop_pair}_vcf2smc.config
echo "dist1=${distidx1}" >> Data/smc++/${pop_pair}/${pop_pair}_vcf2smc.config
echo "dist2=${distidx2}" >> Data/smc++/${pop_pair}/${pop_pair}_vcf2smc.config
echo "outdir=Data/smc++/${pop_pair}" >> Data/smc++/${pop_pair}/${pop_pair}_vcf2smc.config


# -W waits until all jobs finish to continue main script
# conversion step
echo "running vcf2smc"
sbatch -W --array=0-${num_array}%5 Scripts/run_vcf2smc.sh Data/smc++/${pop_pair}/${pop_pair}_vcf2smc.config


# inference step
# WARNING: ASSUMES ESTIMATE HAS BEEN RUN FOR EACH POP
num_pairs=${#pop_pairs[@]}
num_array=$(( num_pairs - 1 ))

mkdir -p Results/smc++/${pop_pair}

splitconfig="Data/smc++/${pop_pair}/${pop_pair}_split.config"
resultsdir="Data/smc++/${pop_pair}"

echo "pop_pair=${pop_pair}" > ${splitconfig}
echo "pop1=${pop1}" >> ${splitconfig}
echo "pop2=${pop2}" >> ${splitconfig}
echo "pop1estimate=Results/smc++/${pop1}_100k_sites_rp5" >> ${splitconfig}
echo "pop2estimate=Results/smc++/${pop2}_100k_sites" >> ${splitconfig}
echo "splitfiles=Data/smc++/${pop_pair}" >> ${splitconfig}
echo "outdir=${resultsdir}" >> ${splitconfig}

echo "running inference"
sbatch -W --partition=dschridelab Scripts/run_split.sh Data/smc++/${pop_pair}/${pop_pair}_split.config
