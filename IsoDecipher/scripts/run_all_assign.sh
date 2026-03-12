#!/bin/bash

EXPS=("exp93" "exp97" "exp100" "exp105" "exp106")
# Use a path relative to where you are running the script
LOG_DIR="results/logs"

# The -p flag creates the directory if it doesn't exist
mkdir -p "$LOG_DIR"

for EXP in "${EXPS[@]}"; do
    echo "Starting $EXP at $(date)"
    # FIX: Ensure there is no "/" before $EXP unless you mean the root
    python IsoDecipher/scripts/assign_reads.py \
        --bam data/${EXP}/possorted_genome_bam.bam \
        --panel results/panel_features.csv \
        --out results/${EXP}_isoform_counts.csv > "${LOG_DIR}/${EXP}.log" 2>&1 & 
done

wait
echo "All samples processed at $(date)."