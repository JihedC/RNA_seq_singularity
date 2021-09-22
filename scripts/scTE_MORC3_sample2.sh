#! /bin/bash

#SBATCH  --job-name=star_solo_sample2
#SBATCH --mail-type=ALL
#SBATCH --mail-user j.chouaref@lumc.nl
#SBATCH -t 6:0:0
#SBATCH --mem=60000

module purge
module load genomics/ngs/samtools/1.11/gcc-8.3.1

scTE -i Morc3_sample2Aligned.sortedByCoord.out.bam -o MORC3_sample2 -x mm10.exclusive.idx -CB CR -UMI UR
scTE -i Morc3_sample2Aligned.sortedByCoord.out.bam -o MORC3_sample2 -x mm10.exclusive.idx --hdf5 True -CB CR -UMI UR
