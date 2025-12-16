#!/bin/bash

#SBATCH --job-name=purc_A
#SBATCH --partition=cpu-g2
#SBATCH --account=tribblelab

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

#SBATCH --array=1-85
#SBATCH --chdir=/mmfs1/gscratch/tribblelab/hemionitidoid_longreads
#SBATCH --export=ALL
#SBATCH -o log/%x_%A_%a.out

module load conda
conda activate purc

# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------
BASE_CONFIG="input/purc_configuration_A.txt"
BATCH_DIR_PARENT="poolafastq"
PYTHON="/mmfs1/gscratch/tribblelab/laufran/purc/purc.py"

# Each array index corresponds to one batch folder.
# Count number of batches automatically:
# sbatch --array=1-N
# ------------------------------------------------------------

# Get batch folder name based on SLURM_ARRAY_TASK_ID
BATCH_NUM=$(printf "%03d" ${SLURM_ARRAY_TASK_ID})
BATCH_FOLDER="${BATCH_DIR_PARENT}/batch_${BATCH_NUM}"
INPUT_FASTQ="${BATCH_FOLDER}/reads.fastq"

# Output naming derived from batch
OUT_PREFIX="poolA_b${BATCH_NUM}"
OUT_DIR="purc_poolA_b${BATCH_NUM}"
CONFIG_OUT="purc_config_b${BATCH_NUM}.conf"

echo "Running PURC for batch ${BATCH_NUM}"
echo "Input FASTQ: ${INPUT_FASTQ}"
echo "Output folder: ${OUT_DIR}"
echo "Output config: ${CONFIG_OUT}"

# ------------------------------------------------------------
# Create modified config file
# ------------------------------------------------------------

# Copy base config
cp "${BASE_CONFIG}" "${CONFIG_OUT}"

# Replace values inside the config file
# (simple in-place edits using sed)
sed -i "s|^Input_sequence_file.*|Input_sequence_file = ${INPUT_FASTQ}|" "${CONFIG_OUT}"
sed -i "s|^Output_prefix.*|Output_prefix = ${OUT_PREFIX}|" "${CONFIG_OUT}"
sed -i "s|^Output_folder.*|Output_folder = ${OUT_DIR}|" "${CONFIG_OUT}"
sed -i "s|^Log_file.*|Log_file = log_${BATCH_NUM}|" "${CONFIG_OUT}"

# ------------------------------------------------------------
# Run PURC
# ------------------------------------------------------------
echo "Launching PURC for batch ${BATCH_NUM}..."
${PYTHON} "${CONFIG_OUT}"

echo "Batch ${BATCH_NUM} complete."

conda deactivate
