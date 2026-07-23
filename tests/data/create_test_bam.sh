#!/usr/bin/env bash

set -euo pipefail

REFERENCE_FASTA='./some_sequence.fasta'
GENOME_LENGTH=1500
SUBSEQUENCES_FASTA='./subsequences.fasta'
SUBSEQ_COUNT=5
SUBSEQ_MIN_LEN=150
SUBSEQ_MAX_LEN=500
SAM_FILE='./some_mapping.sam'
BAM_FILE='./some_mapping.srt.bam'


echo "=== Extracting $SUBSEQ_COUNT random subsequences ==="
# Generate random intervals and extract with seqkit
> "$SUBSEQUENCES_FASTA"  # Empty file

for ((i=1; i<=SUBSEQ_COUNT; i++)); do
    # Random start position (1-based, leaving room for length)
    start=$((RANDOM % (${GENOME_LENGTH} - ${SUBSEQ_MAX_LEN} + 1) + 1))
    length=$((RANDOM % ${SUBSEQ_MAX_LEN} + ${SUBSEQ_MIN_LEN}))  # At least 50bp
    end=$((${start} + ${length} - 1))
    
    # Use seqkit subseq to extract region
    seqkit subseq -r "${start}:${end}" "${REFERENCE_FASTA}"
done | seqkit rename -1 > "${SUBSEQUENCES_FASTA}"


echo "=== Mapping subsequences to original with minimap2 ==="
# Minimap2: -a for SAM output, -x asm5 for assembly-to-assembly
minimap2 -a -x asm5 "${REFERENCE_FASTA}" "${SUBSEQUENCES_FASTA}" > "${SAM_FILE}"

echo "=== Converting to sorted BAM ==="
# Convert SAM to BAM, sort, and index
samtools view -O BAM "${SAM_FILE}" \
    | samtools sort -O BAM -o "${BAM_FILE}"
samtools index "${BAM_FILE}"

rm "${SUBSEQUENCES_FASTA}" ./some_sequence.fasta.seqkit.fai

exit 0
