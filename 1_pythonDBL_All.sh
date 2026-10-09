#!/bin/bash
#SBATCH --job-name=all_python_doublets
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --nodes=1

BASE_DIR="/home/cv1u24/GIST/3_Doublet_Removal"

source /iridisfs/i6software/conda/miniconda-py3/etc/profile.d/conda.sh
conda activate snpatho_env

for SAMPLE_DIR in "${BASE_DIR}"/PATH0*
do
    if [ -d "$SAMPLE_DIR" ]; then
        SAMPLE=$(basename "$SAMPLE_DIR")
        mkdir -p "${SAMPLE_DIR}/SlurmLogs"
        mkdir -p "${SAMPLE_DIR}/Output"
        echo "Processing ${SAMPLE}"
        cd "$SAMPLE_DIR" || exit 1
        python pythonDoublet.py \
> "SlurmLogs/${SAMPLE}_python_doublets.out" \
            2> "SlurmLogs/${SAMPLE}_python_doublets.err"
    fi
done