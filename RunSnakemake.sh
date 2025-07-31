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
module load snakemake OpenBLAS/0.3.17-GCC-11.2.0 ncurses/6.4-GCCcore-12.3.0 SAMtools BWA/0.7.17-foss-2017b picard/2.20.6-Java-1.8.0_144

# Snakemake
snakemake -c12 --rerun-incomplete
