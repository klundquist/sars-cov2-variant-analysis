#!/bin/bash

# Enable alias expansion
shopt -s expand_aliases

source ~/.bash_aliases

# Script: process_sars_cov2_data.sh
# Description: Download FASTQ data, run QC, trim adapters, and generate a combined QC report.
# Author: Karl Lundquist
# Date: 2024-01-31

# Set variables for directories and files
DATA_DIR="data"
QC_DIR="qc_results"
TRIMMED_DIR="trimmed_data"
ADAPTERS_FILE="adapters.fa"
RUN_ID="ERR11584473"
REF_GENOME="SARS-CoV-2-reference.fasta"  # Update this with the actual reference genome file

# Create necessary directories
mkdir -p $DATA_DIR $QC_DIR $TRIMMED_DIR

# Step 1: Check if adapters.fa exists; if not, create it
if [ ! -f $ADAPTERS_FILE ]; then
    echo "Creating adapters.fa file..."
    echo ">Nextera_Adapter_Read1" > $ADAPTERS_FILE
    echo "CTGTCTCTTATACACATCT" >> $ADAPTERS_FILE
    echo ">Nextera_Adapter_Read2" >> $ADAPTERS_FILE
    echo "AGATGTGTATAAGAGACAG" >> $ADAPTERS_FILE
else
    echo "adapters.fa file already exists."
fi

# Step 2: Download FASTQ data using fastq-dump
fastq-dump --split-files --gzip --outdir $DATA_DIR $RUN_ID

# Step 3: Run FastQC on the raw data
fastqc $DATA_DIR/*.fastq.gz -o $QC_DIR

# Step 4: Trim adapters and low-quality bases using Trimmomatic
trimmomatic PE \
    $DATA_DIR/ERR11584473_1.fastq.gz $DATA_DIR/ERR11584473_2.fastq.gz \
    $TRIMMED_DIR/trimmed_forward_paired.fastq.gz $TRIMMED_DIR/trimmed_forward_unpaired.fastq.gz \
    $TRIMMED_DIR/trimmed_reverse_paired.fastq.gz $TRIMMED_DIR/trimmed_reverse_unpaired.fastq.gz \
    ILLUMINACLIP:$ADAPTERS_FILE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

# Step 5: Run FastQC on the trimmed data
fastqc $TRIMMED_DIR/*.fastq.gz -o $QC_DIR

# Step 6: Combine all QC reports using MultiQC
multiqc $QC_DIR -o $QC_DIR

# Step 7: Align trimmed reads to the reference genome using BWA
bwa index $REF_GENOME
bwa mem $REF_GENOME $TRIMMED_DIR/trimmed_forward_paired.fastq.gz $TRIMMED_DIR/trimmed_reverse_paired.fastq.gz > aligned_reads.sam

# Step 8: Convert, sort, and index BAM file
samtools view -Sb aligned_reads.sam > aligned_reads.bam
samtools sort aligned_reads.bam -o sorted_reads.bam
samtools index sorted_reads.bam

# Step 9: Variant calling with FreeBayes
freebayes -f $REF_GENOME sorted_reads.bam > raw_variants.vcf

# Step 10: Variant filtering with GATK
gatk CreateSequenceDictionary -R SARS-CoV-2-reference.fasta -O SARS-CoV-2-reference.dict
gatk VariantFiltration -R $REF_GENOME -V raw_variants.vcf -O filtered_variants.vcf --filter-expression "QD < 2.0" --filter-name "LowQD"

# Step 11: Variant annotation with SnpEff
snpEff ann NC_045512.2 filtered_variants.vcf > annotated_variants.vcf

# Use vcftools to format annotated variants
vcftools --vcf annotated_variants.vcf --get-INFO ANN --out annotations

# Further filter to get a file with one mutation on each line
awk -F '\t' 'NR>1 {split($5, ann, ","); for(i in ann) {split(ann[i], a, "|"); print $1, $2, $3, $4, a[2], a[4], a[10]}}' annotations.INFO > gene_annotations.txt


echo "Analysis complete. Check the $QC_DIR directory for the combined MultiQC report and the output directory for aligned and variant files."
