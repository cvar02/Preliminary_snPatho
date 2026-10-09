#!/bin/bash
#SBATCH --job-name=DoubletAnalysis
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --output=/home/cv1u24/GIST/3_Doublet_Removal/SlurmLogs/DoubletAnalysis_%j.out
#SBATCH --error=/home/cv1u24/GIST/3_Doublet_Removal/SlurmLogs/DoubletAnalysisp_%j.err

source /iridisfs/i6software/conda/miniconda-py3/etc/profile.d/conda.sh
conda activate snpatho_env

cd /home/cv1u24/GIST/3_Doublet_Removal/Scripts || exit 1
python -u 3_DBL_Investigation.py