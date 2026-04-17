#!/usr/bin/env bash

set -euo pipefail


# in_fasta='/mnt/data/Max/genomes/B1622/fasta/v1/Rhodoglobus_sp_B-1622_v1.fa'
in_fasta='/mnt/data/Max/genomes/Decoy_hub/Lactococcus_B-1834_v1.fasta.gz'
workdir='/mnt/data/Max/tmp/minimap_test'

min_len=50
max_len=500
start_coord=2000

if [[ ! -d "${workdir}" ]]; then
    mkdir -pv "${workdir}"
fi


reads_file="${workdir}/reads.fasta"
srt_bam="${workdir}/reads.srt.bam"

for (( len=${min_len}; len<=${max_len}; len++ )); do
    end_coord=$(( ${start_coord} + ${len} - 1 ))
    seq_name="read_len_${len}"
    echo ">${seq_name}"
    seqkit head -n 1 "${in_fasta}" \
        | seqkit subseq --quiet -r "${start_coord}:${end_coord}" \
        | seqkit seq -s
done | seqkit seq -w 100 > "${reads_file}"

echo "${reads_file}"

/home/cager-max/Soft/minimap2-2.30_x64-linux/minimap2 \
    -ax asm5 \
    --eqx \
    -t 1 \
    "${in_fasta}" \
    "${reads_file}" \
| samtools view -O BAM \
| samtools sort -O BAM \
    -o "${srt_bam}"

samtools index "${srt_bam}"

echo "${srt_bam}"


exit 0
