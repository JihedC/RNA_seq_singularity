#!/bin/bash
#SBATCH --job-name=Htseq_count
#SBATCH --time=60:00:00 # hh:mm:ss
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=60000 # 40G
#SBATCH --mail-user=j.chouaref@lumc.nl
#SBATCH --mail-type=ALL
#
#SBATCH --comment=Devepi


htseq-count --format=bam --idattr=transcript_id results/mapped/WT1/WT1.sorted.bam \
  results/mapped/WT2/WT2.sorted.bam \
  results/mapped/WT3/WT3.sorted.bam \
  results/mapped/KO1/KO1.sorted.bam \
  results/mapped/KO2/KO2.sorted.bam \
  results/mapped/KO3/KO3.sorted.bam \
  GRCm38_GENCODE_rmsk_TE.gtf > htseq_Morc3.counts.txt
