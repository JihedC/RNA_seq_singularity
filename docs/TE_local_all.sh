#!/bin/bash
#SBATCH --job-name=TElocal_all
#SBATCH --time=60:00:00 # hh:mm:ss
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=60000 # 40G
#SBATCH --mail-user=j.chouaref@lumc.nl
#SBATCH --mail-type=ALL
#
#SBATCH --comment=Devepi

#Found on github: https://github.com/mhammell-laboratory/TElocal
#TElocal --sortByPos -b RNAseq.bam --GTF gene_annots.gtf --TE te_annots.locInd --project sample_sorted_test

TElocal --sortByPos -b ../results/mapped/WT1/WT1.sorted.bam results/mapped/WT2/WT2.sorted.bam \
results/mapped/WT3/WT3.sorted.bam \
results/mapped/KO1/KO1.sorted.bam \
results/mapped/KO2/KO2.sorted.bam \
results/mapped/KO3/KO3.sorted.bam \
--GTF ../gencode.vM20.annotation.gtf --TE ../mm10_rmsk_TE.gtf.locInd --project test_te_local
