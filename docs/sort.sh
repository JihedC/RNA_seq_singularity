#!/bin/bash
#SBATCH --job-name=sort
#SBATCH --time=60:00:00 # hh:mm:ss
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=60000 # 40G
#SBATCH --mail-user=j.chouaref@lumc.nl
#SBATCH --mail-type=ALL
#
#SBATCH --comment=Devepi

samtools sort -n -o WT_sorted.bam WTAligned.out.bam


samtools sort -n -o KO_sorted.bam KOAligned.out.bam
