#!/bin/bash
#SBATCH --job-name=cellbender_PATH01
#SBATCH --output=/home/cv1u24/GIST/1_Raw_Processing/CellBender/PATH01/SlurmLogs/PATH01.out
#SBATCH --error=/home/cv1u24/GIST/1_Raw_Processing/CellBender/PATH01/SlurmLogs/PATH01.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --nodes=1

#Activate conda environment with python 3.7 (old version, but newest version that is compatible with Cellbender
source /iridisfs/i6software/conda/miniconda-py3/etc/profile.d/conda.sh
conda activate cellbender_py37

#Make needed directories

mkdir -p /home/cv1u24/GIST/1_Raw_Processing/CellBender
mkdir -p /home/cv1u24/GIST/1_Raw_Processing/CellBender/PATH01
mkdir -p /home/cv1u24/GIST/1_Raw_Processing/CellBender/PATH01/SlurmLogs
mkdir -p /home/cv1u24/GIST/1_Raw_Processing/CellBender/PATH01/CellBender_Outputs

#Run cellbender on default settings
cellbender remove-background --input /scratch/cv1u24/GIST/snPatho/PATH01/outs/raw_feature_bc_matrix.h5 --output /home/cv1u24/GIST/1_Raw_Processing/CellBender/PATH01/CellBender_Outputs/PATH01_cellbender.h5

