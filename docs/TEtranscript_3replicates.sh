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


TEtranscripts -t results/mapped/KO1/KO1.sorted.bam \
results/mapped/KO2/KO2.sorted.bam \
results/mapped/KO3/KO3.sorted.bam \
-c results/mapped/WT1/WT1.sorted.bam \
results/mapped/WT2/WT2.sorted.bam \
results/mapped/WT3/WT3.sorted.bam \
--TE GRCm38_Ensembl_rmsk_TE.gtf \
--stranded no \
--format BAM --mode multi \
--GTF gencode.vM20.annotation.gtf  \
--project Morc3_RNA
