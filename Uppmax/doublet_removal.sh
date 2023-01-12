#!/bin/bash -l
#SBATCH -A snic2022-22-723
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:40:00
#SBATCH -J run_python_doublet_removal
#SBATCH --mail-type=ALL
#SBATCH --mail-user Matilda.maxwell.5390@student.uu.se

#modules
module load conda

#commands
source conda_init.sh
conda activate /proj/snic2022-23-388/private/doublet_removal_env

python doublet_removal.py -i /proj/snic2022-23-388/private/data/processed/count_output/$1/outs/filtered_feature_bc_matrix/ -o /proj/snic2022-23-388/private/data/processed/doublet_removal_output/Sample_$1

