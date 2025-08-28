#!/bin/bash
#SBATCH --job-name=align_primer
#SBATCH --output=output_%A_%a.out
#SBATCH --error=error_%A_%a.err
#SBATCH --mem=8G
#SBATCH --array=0-124
#SBATCH --cpus-per-task=8
#SBATCH --partition=compsci

# Parallelize alignment jobs
SEQ_FILES=(Amplicon_1_R2_rc_*)
SEQ_FILE=${SEQ_FILES[$SLURM_ARRAY_TASK_ID]}

# Change to personal conda environment
source /home/users/hiw4/miniconda3/etc/profile.d/conda.sh
conda activate spec_seq_env

python align_spec_seq_primers.py "$SEQ_FILE"
