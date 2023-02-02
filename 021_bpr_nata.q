#!/bin/bash
#
#SBATCH --job-name=bpr_nata
#SBATCH --mail-user=kem81@duke.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=9024
#SBATCH --time=99:00:00
#SBATCH --partition=common
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=slurm-R-job-%J.stdout
#SBATCH --error=slurm-R-job-%J.stderr

module load R/4.1.1

R CMD BATCH 021_bpr_nata.R

module unload R/4.1.1