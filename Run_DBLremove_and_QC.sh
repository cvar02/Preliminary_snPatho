#!/bin/bash
#SBATCH --job-name=DoubletAnalysis
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --output=/home/cv1u24/GIST/3_Doublet_Removal/SlurmLogs/Doublet_Remove_and_QC_%j.out
#SBATCH --error=/home/cv1u24/GIST/3_Doublet_Removal/SlurmLogs/Doublet_Remove_and_QC_%j.err

source /iridisfs/i6software/conda/miniconda-py3/etc/profile.d/conda.sh
conda activate snpatho_env

cd /home/cv1u24/GIST/3_Doublet_Removal/Scripts || exit 1
python -u 4_DBL_Remove_and_QC.py