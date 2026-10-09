#!/bin/bash
#SBATCH --job-name=all_scdblfinder
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --time=48:00:00

BASE_DIR="/home/cv1u24/GIST/3_Doublet_Removal"

source /iridisfs/i6software/conda/miniconda-py3/etc/profile.d/conda.sh
conda activate scdblfinder_env

for SAMPLE_DIR in "${BASE_DIR}"/PATH0*
do
    if [ -d "$SAMPLE_DIR" ]; then
        SAMPLE=$(basename "$SAMPLE_DIR")
        mkdir -p "${SAMPLE_DIR}/SlurmLogs"
        mkdir -p "${SAMPLE_DIR}/Output"
        echo "Processing ${SAMPLE}"
        cd "$SAMPLE_DIR" || exit 1

        Rscript scDBLFinder.R \
> "SlurmLogs/${SAMPLE}_scdblfinder.out" \
            2> "SlurmLogs/${SAMPLE}_scdblfinder.err"
    fi
done
 
