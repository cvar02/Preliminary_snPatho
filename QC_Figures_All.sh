#!/bin/bash
#SBATCH --job-name=QC_all
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --nodes=1

BASE_DIR="/home/cv1u24/GIST/2_QC_Figures"
source /iridisfs/i6software/conda/miniconda-py3/etc/profile.d/conda.sh

conda activate snpatho_env

for SAMPLE_DIR in "$BASE_DIR"/PATH0*; do
   SAMPLE=$(basename "$SAMPLE_DIR")
   LOG_DIR="$SAMPLE_DIR/SlurmLogs"
   mkdir -p "$LOG_DIR"
   echo "Submitting ${SAMPLE}"
   sbatch \
       --job-name="qc_${SAMPLE}" \
       --output="${LOG_DIR}/${SAMPLE}_qc.out" \
       --error="${LOG_DIR}/${SAMPLE}_qc.err" \
       --time=04:00:00 \
       --cpus-per-task=4 \
       --mem=16G \
       --nodes=1 \
       --wrap="cd '$SAMPLE_DIR' && python QCandFigures.py"
done