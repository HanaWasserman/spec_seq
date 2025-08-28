#!/bin/bash
#SBATCH --job-name=align_imperfect
#SBATCH --output=output_%A_%a.out
#SBATCH --error=error_%A_%a.err
#SBATCH --mem=8G
#SBATCH --array=0-9
#SBATCH --cpus-per-task=8
#SBATCH --partition=compsci

SEQ_FILES=(Amplicon_2_R2_imperfect_*)
SEQ_FILE=${SEQ_FILES[$SLURM_ARRAY_TASK_ID]}
LIB_FILE="ETS_library_30_18_1_24_1_18_30_11152024.txt"

source /home/users/hiw4/miniconda3/etc/profile.d/conda.sh

conda activate spec_seq_env

python align_imperfect.py "$SEQ_FILE" "$LIB_FILE"
