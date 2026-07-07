#!/bin/bash
#SBATCH --job-name=ccs_conversion
#SBATCH --partition=cpu-g2
#SBATCH --account=tribblelab
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --chdir=/mmfs1/gscratch/tribblelab/hemionitidoid_longreads
#SBATCH --export=ALL
#SBATCH -o log/ccs_%j.out
#SBATCH -e log/ccs_%j.err

# ------------------------------------------------------------
# User settings - UPDATE THESE PATHS
# ------------------------------------------------------------
# Path to your demultiplexed subreads bam files
SUBREADS_POOLA="/gscratch/tribblelab/hemionitidoid_longreads/input/m64044_221017_183110.subreads.demux.pentpoolA--pentpoolA.bam"
SUBREADS_POOLB="/gscratch/tribblelab/hemionitidoid_longreads/input/m64044_221017_183110.subreads.demux.pentpoolB--pentpoolB.bam"

# Output directory
OUTPUT_DIR="/mmfs1/gscratch/tribblelab/hemionitidoid_longreads/input"

# Number of threads (should match --cpus-per-task)
THREADS=16

module load conda
conda activate pacbio

# ------------------------------------------------------------
# Create output directory
# ------------------------------------------------------------
mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

# ------------------------------------------------------------
# Run CCS on poolA
# ------------------------------------------------------------
echo "Starting CCS for poolA at $(date)"
ccs ${SUBREADS_POOLA} \
    poolA.ccs.bam \
    --min-rq 0.9 \
    --min-passes 3 \
    --num-threads ${THREADS} \
    --log-level INFO \
    --log-file poolA.ccs.log

echo "CCS for poolA completed at $(date)"

# Convert poolA to fastq
echo "Converting poolA to fastq..."
bam2fastq -o poolA.ccs poolA.ccs.bam

# This creates poolA.ccs.fastq.gz - decompress it
#gunzip poolA.ccs.fastq.gz

# ------------------------------------------------------------
# Run CCS on poolB
# ------------------------------------------------------------
echo "Starting CCS for poolB at $(date)"
ccs ${SUBREADS_POOLB} \
  poolB.ccs.bam \
    --min-rq 0.9 \
    --min-passes 3 \
    --num-threads ${THREADS} \
    --log-level INFO \
    --log-file poolB.ccs.log

echo "CCS for poolB completed at $(date)"

# Convert poolB to fastq
echo "Converting poolB to fastq..."
bam2fastq -o poolB.ccs poolB.ccs.bam

# This creates poolB.ccs.fastq.gz - decompress it
gunzip poolB.ccs.fastq.gz

conda deactivate