#!/bin/bash
#SBATCH --output=Logs/Sample1-Sample2.out
#SBATCH --partition={Node}
#SBATCH --mem=500G
#SBATCH -c 12

# Make directories
cd /Alice/Snakemake-WGS/
if [ ! -d ./Samtools-Output/ ]; then mkdir -p ./Samtools-Output/; fi
if [ ! -d ./VCF-Output/ ]; then mkdir -p ./VCF-Output/; fi
if [ ! -d ./VEP-Output/ ]; then mkdir -p ./VEP-Output/; fi

# Load modules
module load snakemake OpenBLAS ncurses SAMtools BWA picard

# Snakemake
snakemake -c12 --rerun-incomplete
