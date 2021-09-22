#! /bin/bash

#SBATCH  --job-name=star_solo
#SBATCH --mail-type=ALL
#SBATCH --mail-user j.chouaref@lumc.nl
#SBATCH -t 1:0:0
#SBATCH --mem=60000

READS="/exports/archive/hg-groep-vandermaarel/Lucia/single_cell_projects/Thymus_scRNAseq/Sample_2"
INDEX="/exports/humgen/jihed/TE_Transcript/temp/genome"

module purge
module load genomics/ngs/samtools/1.11/gcc-8.3.1

STAR --genomeDir $INDEX \
    --readFilesIn $READS/Sample_2_S2_L003_R2_001.fastq.gz $READS/Sample_2_S2_L003_R1_001.fastq.gz \
    --soloUMIlen 12 --soloType CB_UMI_Simple --soloCBwhitelist 3M-february-2018.txt  \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI nM AS CR UR CB UB \
    --readFilesCommand zcat \
    --outFilterMultimapNmax 100 \
    --winAnchorMultimapNmax 100 \
    --outMultimapperOrder Random \
    --runRNGseed 777 --outSAMmultNmax 1\
    --outFileNamePrefix Morc3_sample2
