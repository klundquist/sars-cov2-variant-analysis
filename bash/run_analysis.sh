#!/bin/bash

# Script: run_analysis.sh
# Description: Wrapper script for SARS-CoV-2 variant analysis pipeline
# Author: Karl Lundquist
# Date: 2025-03-20

# Default values
DATA_DIR="../data"
OUTPUT_DIR="../results"
RUN_ID="ERR11584473"
REF_GENOME="${DATA_DIR}/SARS-CoV-2-reference.fasta"
ADAPTERS_FILE="${DATA_DIR}/adapters.fa"

# Help function
function show_help {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -h, --help                  Show this help message"
    echo "  -d, --data-dir DIR          Input data directory (default: ../data)"
    echo "  -o, --output-dir DIR        Output directory (default: ../results)"
    echo "  -r, --run-id ID             SRA run ID (default: ERR11584473)"
    echo "  -g, --reference FILE        Reference genome file (default: ../data/SARS-CoV-2-reference.fasta)"
    echo "  -a, --adapters FILE         Adapters file (default: ../data/adapters.fa)"
    exit 0
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -h|--help)
            show_help
            ;;
        -d|--data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -r|--run-id)
            RUN_ID="$2"
            shift 2
            ;;
        -g|--reference)
            REF_GENOME="$2"
            shift 2
            ;;
        -a|--adapters)
            ADAPTERS_FILE="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            ;;
    esac
done

# Create output directories
QC_DIR="${OUTPUT_DIR}/qc_results"
TRIMMED_DIR="${OUTPUT_DIR}/trimmed_data"
mkdir -p $DATA_DIR $QC_DIR $TRIMMED_DIR

# Export variables for the pipeline script
export DATA_DIR
export QC_DIR
export TRIMMED_DIR
export ADAPTERS_FILE
export RUN_ID
export REF_GENOME
export OUTPUT_DIR

# Check if required files exist
if [ ! -f "$ADAPTERS_FILE" ]; then
    echo "Adapters file not found. Creating default adapters.fa."
    echo ">Nextera_Adapter_Read1" > "$ADAPTERS_FILE"
    echo "CTGTCTCTTATACACATCT" >> "$ADAPTERS_FILE"
    echo ">Nextera_Adapter_Read2" >> "$ADAPTERS_FILE"
    echo "AGATGTGTATAAGAGACAG" >> "$ADAPTERS_FILE"
fi

# Execute the pipeline script
echo "Starting SARS-CoV-2 variant analysis pipeline..."
SCRIPT_DIR=$(dirname "$0")
cd "$OUTPUT_DIR" && bash "$SCRIPT_DIR/pipeline.sh"

echo "Analysis complete. Results are available in $OUTPUT_DIR"