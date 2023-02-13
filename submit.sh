#!/bin/sh
#SBATCH -t 01:00:00
#SBATCH -A mulif007a
#SBATCH -N1
#SBATCH -p ProdQ
#SBATCH --job-name=analysis


cd $SLURM_SUBMIT_DIR

module load conda/2

source activate /ichec/work/mulif007a/ojas/envs/analysis

python prepare.py
