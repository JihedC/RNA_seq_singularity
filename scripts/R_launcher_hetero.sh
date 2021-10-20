#! /bin/bash

#SBATCH  --job-name=Seurat_hetero
#SBATCH --mail-type=ALL
#SBATCH --mail-user j.chouaref@lumc.nl
#SBATCH -t 2:0:0
#SBATCH --mem=100000

conda activate seurat 
Rscript Seurat_heterozygous.R
