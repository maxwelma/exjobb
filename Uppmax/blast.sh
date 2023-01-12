#!/bin/bash -l
#SBATCH -A snic2022-22-723
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 08:00:00
#SBATCH -J blastall
#SBATCH --mail-type=ALL
#SBATCH --mail-user Matilda.maxwell.5390@student.uu.se
# Load modules
module load bioinfo-tools
module load blast

cat   */*.fastq > all_seqs.fastq

sed -n '1~4s/^@/>/p;2~4p' all_seqs.fastq > all_seqs.fasta

makeblastdb -in all_seqs.fasta -dbtype nucl -out all_seqs/all_seqs

blastn -task blastn-short -db all_seqs/all_seqs  -query query.fasta -out outfileblast

