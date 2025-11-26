#!/bin/bash

#SBATCH --job-name=purc_A
#SBATCH --partition=cpu-g2
#SBATCH --account=tribblelab

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=480G
#SBATCH --time=48:00:00

#SBATCH --chdir=/mmfs1/gscratch/tribblelab/hemionitidoid_longreads
#SBATCH --export=/mmfs1/gscratch/tribblelab/hemionitidoid_longreads
#SBATCH -o log/%x_%j.out

conda init
conda activate purc

/mmfs1/gscratch/tribblelab/laufran/purc/purc.py scripts/purc_configuration_A.txt