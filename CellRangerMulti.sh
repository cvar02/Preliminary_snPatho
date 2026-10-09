#!/bin/bash
 
#SBATCH --job-name=cellranger_pipeline
#SBATCH --output=/home/cv1u24/GIST/CellRanger/SlurmLogs/%x_%j.out
#SBATCH --error=/home/cv1u24/GIST/CellRanger/SlurmLogs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=48:00:00
 
# -----------------------------
# Paths and settings
# -----------------------------
OUTPUT_BASE="/home/cv1u24/GIST/CellRanger/Output"
CONFIG_BASE="/home/cv1u24/GIST/References/snPatho"
SAMPLES=("PATH01" "PATH02" "PATH03" "PATH04")
 
mkdir -p "$OUTPUT_BASE"
mkdir -p "/home/cv1u24/GIST/CellRanger/SlurmLogs"
 
# -----------------------------
# Load Cell Ranger
# -----------------------------
export PATH=$HOME/GIST/References/cellranger-10.0.0/bin:$PATH
 
# -----------------------------
# Check config CSVs exist
# -----------------------------
echo "Checking config CSVs for all samples..."
 
for SAMPLE in "${SAMPLES[@]}"; do
    NUM="${SAMPLE: -2}"  # extracts 01, 02, 03, 04
    CSV="$CONFIG_BASE/MultiCellRangerPATH${NUM}.csv"
 
    if [ ! -f "$CSV" ]; then
        echo "ERROR: Config CSV not found: $CSV"
        exit 1
    fi
    echo "Found: $CSV"
done
 
echo "All config CSVs found"
 
# -----------------------------
# Run cellranger multi
# -----------------------------
cd "$OUTPUT_BASE"
 
for SAMPLE in "${SAMPLES[@]}"; do
    NUM="${SAMPLE: -2}"
    CSV="$CONFIG_BASE/MultiCellRangerPATH${NUM}.csv"
 
    echo "==============================="
    echo "Running Cell Ranger multi for $SAMPLE"
    echo "==============================="
 
    cellranger multi \
        --id="$SAMPLE" \
        --csv="$CSV" \
        --localcores=40 \
        --localmem=128
done
 
# -----------------------------
# Create aggregation CSV
# -----------------------------
AGGR_CSV="$OUTPUT_BASE/aggregation.csv"
echo "Creating aggregation CSV at $AGGR_CSV..."
echo "sample_id,molecule_h5" > "$AGGR_CSV"
 
for SAMPLE in "${SAMPLES[@]}"; do
    H5_FILE="$OUTPUT_BASE/$SAMPLE/outs/per_sample_outs/$SAMPLE/count/molecule_info.h5"
 
    if [ -f "$H5_FILE" ]; then
        echo "$SAMPLE,$H5_FILE" >> "$AGGR_CSV"
    else
        echo "WARNING: $H5_FILE not found, skipping $SAMPLE"
    fi
done
 
# -----------------------------
# Run aggregation
# -----------------------------
echo "Running Cell Ranger aggregation..."
 
cellranger aggr \
    --id="GIST_AllSamples" \
    --csv="$AGGR_CSV" \
    --normalize=mapped
 
echo "Pipeline completed!"
