#!/bin/bash

#SBATCH --job-name=gz_decompress
#SBATCH --partition=cpu-g2
#SBATCH --account=tribblelab

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=480G
#SBATCH --time=6:00:00

#SBATCH --chdir=/mmfs1/gscratch/tribblelab/hemionitidoid_longreads
#SBATCH --export=/mmfs1/gscratch/tribblelab/hemionitidoid_longreads
#SBATCH -o log/%x_%j.out


gunzip /mmfs1/gscratch/tribblelab/hemionitidoid_longreads/input/pentpoolA.fastq.gz
#gunzip /mmfs1/gscratch/tribblelab/hemionitidoid_longreads/input/pentpoolB.fastq.gz
