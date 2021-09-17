#! /bin/bash

#SBATCH  --job-name=star_solo
#SBATCH --mail-type=ALL
#SBATCH --mail-user j.chouaref@lumc.nl
#SBATCH -t 24:00:00
#SBATCH --mem=60000

#/path/to/STAR --genomeDir /path/to/genome/dir/ --readFilesIn ...  [...other parameters...] --soloType ... --soloCBwhitelist ...

READS="/exports/humgen/jihed/scTE/Kernfeld_e12"
INDEX="/exports/humgen/jihed/TE_Transcript/temp/genome"
module purge
module load genomics/ngs/samtools/1.11/gcc-8.3.1

STAR --genomeDir $INDEX \
    --readFilesIn $READS/SRX3461253_T1_2.fastq.gz $READS/SRX3461253_T1_1.fastq.gz \
    --soloUMIlen 12 --soloType CB_UMI_Simple --soloCBwhitelist 3M-february-2018.txt  \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI nM AS CR UR CB UB \
    --readFilesCommand zcat \
    --outFilterMultimapNmax 100 \
    --winAnchorMultimapNmax 100 \
    --outMultimapperOrder Random \
    --runRNGseed 777 --outSAMmultNmax 1\
    --outFileNamePrefix Kernfeld_e12/e12_rep2

cd Kernfeld_e12/
scTE -i e12_rep2Aligned.sortedByCoord.out.bam -o Kernfeld_e12_rep2 -x ../mm10.exclusive.idx --hdf5 True -CB CR -UMI UR
scTE -i e12_rep2Aligned.sortedByCoord.out.bam -o Kernfeld_e12_rep2 -x ../mm10.exclusive.idx -CB CR -UMI UR    