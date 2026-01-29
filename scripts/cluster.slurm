#!/bin/bash

#SBATCH --job-name=purc_recluster
#SBATCH --partition=cpu-g2
#SBATCH --account=tribblelab

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G                     # Increased memory for heavy chimeric checking
#SBATCH --time=48:00:00

#SBATCH --array=0-8
#SBATCH --chdir=/mmfs1/gscratch/tribblelab/hemionitidoid_longreads
#SBATCH --export=ALL
#SBATCH -o log/%x_%A_%a.out

module load conda
conda activate purc

PROJECT_ROOT="/mmfs1/gscratch/tribblelab/hemionitidoid_longreads"
BATCH_DIR_PARENT="${PROJECT_ROOT}/output/dummyrun/purc_poolA_b001"
PYTHON_SCRIPT="${PROJECT_ROOT}/scripts/purc_recluster.py"

#locus list (fixit: need to add APPEFP eventually)
LOCI=(
  CRY2
  GAPCPSH
  IBR3
  INFA
  MATK
  PGIC
  PSBM
  RPL2
  SQD1
)

LOCUS="${LOCI[$SLURM_ARRAY_TASK_ID]}"

echo "Reclustering locus: ${LOCUS}"

python ${PYTHON_SCRIPT} \
  -f "${BATCH_DIR_PARENT}/${LOCUS}/${LOCUS}.fa" \
  -o "${BATCH_DIR_PARENT}/${LOCUS}/${LOCUS}_recluster" \
  -c 0.997 0.995 0.99 0.997 \
  -s 1 4 \
  --clean \
  -m OTU

conda deactivate