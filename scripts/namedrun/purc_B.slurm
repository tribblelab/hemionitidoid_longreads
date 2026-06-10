#!/bin/bash
#SBATCH --job-name=purc_poolB
#SBATCH --partition=cpu-g2
#SBATCH --account=tribblelab
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=38G
#SBATCH --time=72:00:00
#SBATCH --chdir=/mmfs1/gscratch/tribblelab/hemionitidoid_longreads
#SBATCH --export=ALL
#SBATCH -o log/%x_%j.out
#SBATCH -e log/%x_%j.err

module load conda
conda activate purc

# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------
BASE_CONFIG="/mmfs1/gscratch/tribblelab/hemionitidoid_longreads/input/namedrefsrun/purc_configuration_B.txt"
PROJECT_ROOT="/mmfs1/gscratch/tribblelab/hemionitidoid_longreads"
OUTPUT_ROOT="${PROJECT_ROOT}/output"
PYTHON="${PROJECT_ROOT}/../laufran/purc/purc.py"

# ------------------------------------------------------------
# CCS input file
# ------------------------------------------------------------
INPUT_FASTQ="${PROJECT_ROOT}/input/poolB.ccs.fastq"
OUT_PREFIX="poolB"
OUT_DIR="${OUTPUT_ROOT}/purc_poolB"
CONFIG_OUT="${PROJECT_ROOT}/input/purc_config_poolB.conf"

echo "Running PURC on CCS reads"
echo "Input FASTQ: ${INPUT_FASTQ}"
echo "Output folder: ${OUT_DIR}"
echo "Output config: ${CONFIG_OUT}"

# ------------------------------------------------------------
# Create modified config file
# ------------------------------------------------------------
# Copy base config
cp "${BASE_CONFIG}" "${CONFIG_OUT}"

# Replace values inside the config file
sed -i "s|^Input_sequence_file.*|Input_sequence_file = ${INPUT_FASTQ}|" "${CONFIG_OUT}"
sed -i "s|^Output_prefix.*|Output_prefix = ${OUT_PREFIX}|" "${CONFIG_OUT}"
sed -i "s|^Output_folder.*|Output_folder = ${OUT_DIR}|" "${CONFIG_OUT}"
sed -i "s|^Log_file.*|Log_file = log_poolA_ccs|" "${CONFIG_OUT}"

# ------------------------------------------------------------
# Run PURC
# ------------------------------------------------------------
echo "Launching PURC on CCS reads at $(date)..."
${PYTHON} "${CONFIG_OUT}"

echo "PURC complete at $(date)"
echo "Output directory: ${OUT_DIR}"

conda deactivate