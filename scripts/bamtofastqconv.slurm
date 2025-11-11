#!/bin/bash

#SBATCH --job-name=fastq_conv
#SBATCH --partition=cpu-g2
#SBATCH --account=tribblelab

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=480G
#SBATCH --time=6:00:00

#SBATCH --chdir=/mmfs1/gscratch/tribblelab/hemionitidoid_longreads
#SBATCH --export=/mmfs1/gscratch/tribblelab/hemionitidoid_longreads
#SBATCH -o log/%x_%j.out

conda init
conda activate pacbio

/mmfs1/gscratch/tribblelab/laufran/miniconda3/envs/pacbio/bin/bam2fastq -o input/pentpoolA input/m64044_221017_183110.subreads.demux.pentpoolA--pentpoolA.bam
/mmfs1/gscratch/tribblelab/laufran/miniconda3/envs/pacbio/bin/bam2fastq -o input/pentpoolB input/m64044_221017_183110.subreads.demux.pentpoolB--pentpoolB.bam
