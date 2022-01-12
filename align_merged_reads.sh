#!/bin/bash

source .env
ANALYSIS_DIR="$1"
REFERENCE_SEQ="${PROJECTPATH}/bwa_reference/indexed/reference_pseudoORF.fa"

# index reference genome
${BWA_PATH} index "${REFERENCE_SEQ}"

# aligns reads
${BWA_PATH} mem "${REFERENCE_SEQ}" "${ANALYSIS_DIR}/merged_reads.fastq" | \
samtools sort -@8 -o "${ANALYSIS_DIR}/aligned_reads_sorted.bam" -

# lists allelic frequencies
bcftools mpileup -f "${REFERENCE_SEQ}" \
  -a FORMAT/AD,INFO/AD \
  "${ANALYSIS_DIR}/aligned_reads_sorted.bam" > \
  "${ANALYSIS_DIR}/var.vcf"
