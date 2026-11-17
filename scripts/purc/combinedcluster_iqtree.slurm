#!/bin/bash
#SBATCH --job-name=pipeline_submit
#SBATCH --partition=cpu-g2
#SBATCH --account=tribblelab
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --chdir=/mmfs1/gscratch/tribblelab/hemionitidoid_longreads
#SBATCH -o log/pipeline_submit_%j.out

# =============================================================================
# Pipeline launcher — submits cluster job, then chains IQ-TREE with dependency
# Submit with: sbatch cluster_iqtree.slurm
# =============================================================================

POOL="poolB"       # ← only change this
RUN_NAME="june3"  # ← and this

# -----------------------------------------------------------------------------
# Step 1: Submit clustering
# -----------------------------------------------------------------------------
CLUSTER_JID=$(sbatch --parsable --export=ALL,POOL=$POOL,RUN_NAME=$RUN_NAME <<'CLUSTER'
#!/bin/bash
#SBATCH --job-name=recluster
#SBATCH --partition=cpu-g2
#SBATCH --account=tribblelab
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=48:00:00
#SBATCH --array=0-8
#SBATCH --chdir=/mmfs1/gscratch/tribblelab/hemionitidoid_longreads
#SBATCH -o log/recluster_%A_%a.out

module load conda
conda activate purc

PROJECT_ROOT="/mmfs1/gscratch/tribblelab/hemionitidoid_longreads"
BATCH_DIR_PARENT="${PROJECT_ROOT}/output/${RUN_NAME}/purc_${POOL}"
PYTHON_SCRIPT="${PROJECT_ROOT}/scripts/purc_recluster.py"

LOCI=(APPEFP CRY2 IBR3 INFA MATK PGIC PSBM RPL2 SQD1)
LOCUS="${LOCI[$SLURM_ARRAY_TASK_ID]}"

echo "Pool: ${POOL}  Locus: ${LOCUS}"
echo "Input dir: ${BATCH_DIR_PARENT}"

cd "${BATCH_DIR_PARENT}" || { echo "ERROR: could not cd to ${BATCH_DIR_PARENT}"; exit 1; }

python ${PYTHON_SCRIPT} \
  -f "${LOCUS}/${LOCUS}.fa" \
  -o "${LOCUS}/${LOCUS}_recluster" \
  -c 0.997 0.995 0.99 0.997 \
  -s 1 4 \
  --clean \
  -m OTU

conda deactivate
CLUSTER
)

echo "Submitted clustering job: ${CLUSTER_JID}"

# -----------------------------------------------------------------------------
# Step 2: Submit IQ-TREE with dependency on clustering
# -----------------------------------------------------------------------------
IQTREE_JID=$(sbatch --parsable \
    --dependency=afterok:${CLUSTER_JID} \
    --kill-on-invalid-dep=yes \
    --export=ALL,POOL=$POOL,RUN_NAME=$RUN_NAME <<'IQTREE'
#!/bin/bash
#SBATCH --job-name=iqtree_genetrees
#SBATCH --partition=cpu-g2
#SBATCH --account=tribblelab
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --array=0-8
#SBATCH --chdir=/mmfs1/gscratch/tribblelab/hemionitidoid_longreads
#SBATCH --output=log/iqtree_%A_%a.out
#SBATCH --error=log/iqtree_%A_%a.err

module load coenv/iqtree2/2.3.6

BASE="/gscratch/tribblelab/hemionitidoid_longreads/output"
POOL_DIR="${BASE}/${RUN_NAME}/purc_${POOL}"
OUTPUT_BASE="${BASE}/${RUN_NAME}/${POOL}"

GENES=(APPEFP CRY2 IBR3 INFA MATK PGIC PSBM RPL2 SQD1)
GENE="${GENES[$SLURM_ARRAY_TASK_ID]}"

ALIGNMENT="${POOL_DIR}/${GENE}/${GENE}_recluster/${GENE}_OTUs.aligned.fa"

if [[ ! -f "$ALIGNMENT" ]]; then
    echo "ERROR: Alignment not found: ${ALIGNMENT}"
    exit 1
fi

OUTDIR="${OUTPUT_BASE}/${GENE}"
mkdir -p "$OUTDIR"

echo "========================================"
echo "Task ID   : ${SLURM_ARRAY_TASK_ID}"
echo "Pool      : ${POOL}"
echo "Gene      : ${GENE}"
echo "Alignment : ${ALIGNMENT}"
echo "Output    : ${OUTDIR}"
echo "Start     : $(date)"
echo "========================================"

iqtree2 \
    -s "$ALIGNMENT" \
    -m TEST \
    -B 1000 \
    -T "$SLURM_CPUS_PER_TASK" \
    --prefix "${OUTDIR}/${GENE}"

EXIT_CODE=$?
echo "End time : $(date)"
echo "Exit code: ${EXIT_CODE}"
exit $EXIT_CODE
IQTREE
)

echo "Submitted IQ-TREE job: ${IQTREE_JID} (depends on ${CLUSTER_JID})"
echo "Monitor with: squeue -u $USER"
