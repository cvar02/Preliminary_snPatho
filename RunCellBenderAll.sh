#!/bin/bash
#SBATCH --job-name=cellbender_all
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --nodes=1

BASE_DIR="/home/cv1u24/GIST/1_Raw_Processing/CellBender"
for script in "$BASE_DIR"/PATH0*/RunCellBender.sh; do
   if [ -f "$script" ]; then
       sample_dir=$(dirname "$script")
       sample_name=$(basename "$sample_dir")
       echo "Submitting CellBender job for $sample_name"
       sbatch "$script"
   else
       echo "No RunCellBender.sh scripts found."
   fi
done