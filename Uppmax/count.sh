#!/bin/bash -l
#SBATCH -A snic2022-22-723
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 08:00:00
#SBATCH -J cellranger_count
#SBATCH --mail-type=ALL
#SBATCH --mail-user Matilda.maxwell.5390@student.uu.se
# Load modules
module load bioinfo-tools
module load cellranger
#commands
sample=$1
#cellranger count --id=$sample --fastqs=/proj/snic2022-23-388/private/data/unzipped/Sample_$sample --sample=$sample --transcriptome=/proj/snic2022-23-388/private/data/cellranger_count_test/GRCg7b_test_22-09-12
cellranger count --id=$sample --fastqs=/proj/snic2022-23-388/private/data/unzipped/Sample_$sample --sample=$sample --transcriptome=/crex/proj/snic2022-23-388/private/scripts/evas_gen_test 

