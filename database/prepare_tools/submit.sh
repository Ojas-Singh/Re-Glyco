#!/bin/sh
#SBATCH -t 01:00:00
#SBATCH -A mulif007a
#SBATCH -N1
#SBATCH -p DevQ
#SBATCH --job-name=analysis
#SBATCH --mail-user=ojas.singh.2023@mumail.ie
#SBATCH --mail-type=END

cd $SLURM_SUBMIT_DIR

module load conda/2

source activate /ichec/work/mulif007a/ojas/envs/analysis

python main.py
