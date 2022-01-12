#!/bin/bash

source .env
echo "$BWA_DIR/bwa"
TEST_READS_DIR="${DATADIR}/out/1_1pre_GEM_325406w_BE2"
REFERENCE_SEQ="${PROJECTPATH}/bwa_reference/indexed/reference_pseudoORF.fa"

# index reference genome
${BWA_PATH} index "${REFERENCE_SEQ}"

# aligns reads
${BWA_PATH} mem "${REFERENCE_SEQ}" "${TEST_READS_DIR}/merged_reads.fastq" | \
samtools sort -@8 -o "${TEST_READS_DIR}/aligned_reads_sorted.bam" -

# lists allelic frequencies
bcftools mpileup -f "${REFERENCE_SEQ}" \
  -a FORMAT/AD,INFO/AD \
  "${TEST_READS_DIR}/aligned_reads_sorted.bam" > \
  "${TEST_READS_DIR}/var.vcf"
