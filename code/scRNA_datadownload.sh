#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 6-00:00 # Runtime in D-HH:MM
#SBATCH -J scRNA
#SBATCH --output=/ix/djishnu/Trirupa/ABomics.Prj/SigmiR/scripts/scRNAseq_datadownload.out
#SBATCH --mail-type=END,FAIL
#SBATCH --cluster=htc
#SBATCH --mail-user=trc84@pitt.edu
#SBATCH --account=djishnu
#SBATCH --mem=100g


module load gcc/12.2.0
module load r/4.3.0

Rscript /ix/djishnu/Trirupa/ABomics.Prj/scripts/scRNA_datadownload.R
echo done
