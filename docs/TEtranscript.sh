#!/bin/bash
#SBATCH --job-name=TEtranscript
#SBATCH --time=60:00:00 # hh:mm:ss
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=60000 # 40G
#SBATCH --mail-user=j.chouaref@lumc.nl
#SBATCH --mail-type=ALL
#
#SBATCH --comment=Devepi

TEtranscripts -t mapped/WT_sorted.bam -c mapped/KO_sorted.bam --TE GRCm38_Ensembl_rmsk_TE.gtf --stranded no --format BAM --mode multi --GTF gencode.vM20.annotation.gtf --project mytestProject
