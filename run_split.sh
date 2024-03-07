#!/bin/bash
#SBATCH -D /work/users/t/v/tvkent/aedes/
#SBATCH --time=50:00:00
#SBATCH --mem=100G
#SBATCH -J split
#SBATCH -o /work/users/t/v/tvkent/aedes/logs/split-out-%A.%a.txt
#SBATCH -e /work/users/t/v/tvkent/aedes/logs/split-error-%A.%a.txt

set -e
set -u

workdir=/work/users/t/v/tvkent/aedes

# source arguments from config
source $1

#pop_pairs=("SoCal-Rio_Claro" "NorCal-Rio_Claro" "Florida-Rio_Claro" "Exeter-Rio_Claro" "Clovis-Rio_Claro" "SoCal-Senegal" "SoCal-Gabon")
#pop_pairs=("Gabon-USA")
#pop_pair=${pop_pairs[SLURM_ARRAY_TASK_ID]}

pop1=${pop_pair%-*}
pop2=${pop_pair#*-}

module load singularity
#echo ${splitfiles}
#outdir="Data/smc++/${pop_pair}"
#pop1estimate="Results/smc++/${pop1}_100k_sites"
#pop2estimate="Results/smc++/${pop2}_100k_sites"
#splitfiles="Data/smc++/${pop_pair}"

singularity exec --bind /work/users/t/v/tvkent/aedes/ ./smc++ smc++ split -o ${outdir} ${pop1estimate}/model.final.json ${pop2estimate}/model.final.json ${splitfiles}/*smc.gz

singularity exec --bind /work/users/t/v/tvkent/aedes/ ./smc++ smc++ plot ${outdir}/${pop_pair}.pdf ${outdir}/model.final.json
