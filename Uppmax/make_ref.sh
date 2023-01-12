#!/bin/bash -l
#SBATCH -A snic2022-22-723
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 08:00:00
#SBATCH -J cellranger_mkref
#SBATCH --mail-type=ALL
#SBATCH --mail-user Matilda.maxwell.5390@student.uu.se
# Load modules
module load bioinfo-tools
module load cellranger
#commands
#cellranger mkref --genome=GRCg7b_test_22-09-12 --fasta=fasta/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna --genes=genes/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gtf
cellranger mkref --genome=evas_gen_test --fasta=/proj/snic2022-23-388/private/data/genome/ncbi-genomes-2022-09-09/fasta/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna --genes=/proj/snic2022-23-388/private/data/annotation/ncbi-genomes-2022-09-09/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gtf

