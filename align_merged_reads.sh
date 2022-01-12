#!/bin/bash

source .env
echo "$BWA_DIR/bwa"
${BWA_DIR}/bwa index bwa_reference/reference_pseudoORF.fa