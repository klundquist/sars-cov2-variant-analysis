#!/bin/bash

# Script: run_snakemake.sh
# Description: Wrapper script for SARS-CoV-2 variant analysis with Snakemake
# Author: Karl Lundquist
# Date: 2025-03-20

# Default values
CORES=4
DRY_RUN=false
REASON=false

# Help function
function show_help {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -h, --help                  Show this help message"
    echo "  -c, --cores N               Number of cores to use (default: 4)"
    echo "  -n, --dry-run               Perform a dry run (don't execute any commands)"
    echo "  -r, --reason                Print the reason for each executed rule"
    exit 0
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -h|--help)
            show_help
            ;;
        -c|--cores)
            CORES="$2"
            shift 2
            ;;
        -n|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -r|--reason)
            REASON=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            ;;
    esac
done

# Build command
CMD="snakemake --cores $CORES"

if [ "$DRY_RUN" = true ]; then
    CMD="$CMD -n"
fi

if [ "$REASON" = true ]; then
    CMD="$CMD -r"
fi

# Execute Snakemake
echo "Starting SARS-CoV-2 variant analysis pipeline with Snakemake..."
echo "Command: $CMD"
eval $CMD

exit_code=$?
if [ $exit_code -eq 0 ]; then
    echo "Snakemake analysis complete."
else
    echo "Snakemake analysis failed with exit code $exit_code."
fi