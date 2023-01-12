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
module load bwa
module load blast

while getopts q:i:o: flag
do
    case "${flag}" in
        q) query=${OPTARG};;
        i) infile=${OPTARG};;
        o) outfile=${OPTARG};;
    esac
done

#bwa mem -w 1  $query $infile > "$outfile"

sed -n '1~4s/^@/>/p;2~4p' $infile > $infile.fasta

makeblastdb -in $infile.fasta -dbtype nucl -out blast_dbs/$infile
blastn -task blastn-short -db blast_dbs/$infile  -query $query > blast_results/$outfile

