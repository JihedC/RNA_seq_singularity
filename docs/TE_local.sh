#!/bin/bash
#SBATCH --job-name=TElocal
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

TElocal --sortByPos -b ../results/mapped/WT1/WT1.sorted.bam --GTF ../gencode.vM20.annotation.gtf --TE ../ --project test_te_local
