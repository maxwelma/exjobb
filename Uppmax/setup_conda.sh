#!/bin/bash -l
#SBATCH -A snic2022-22-723
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 04:00:00
#SBATCH -J create_conda_env
#SBATCH --mail-type=ALL
#SBATCH --mail-user Matilda.maxwell.5390@student.uu.se
module load conda

echo "Doublet removal python3 conda enviroment doublet_removal_env"
echo "Setting up..."

export CONDA_ENVS_PATH=/proj/snic2022-23-388/private/

###create the conda environment 
conda env create -n doublet_removal_env --file doublet_removal_env.yml

