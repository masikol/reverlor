#!/usr/bin/env bash

set -euo pipefail

FINDER='repeatscout' # choice
# FINDER='minimap2' # choice

MASK_ALL_BUT_REPEATS=0 # choice
# MASK_ALL_BUT_REPEATS=1 # choice

N_THREADS=8

WORKDIR='/mnt/data/Max/repos/reverlor'
reverlor_find_fpath="${WORKDIR}/reverlor/reverlor_find.py"
genome_fasta="${WORKDIR}/param_selection/data/Mycoplasma_mycoides_JCVI-syn3.0.fasta"
pipeline_dir="${WORKDIR}/param_selection/test_repeat_detection"

pipeline_workdir="${pipeline_dir}/workdir_${FINDER}"
find_repeats_out_root="${pipeline_workdir}/find_repeats_results"
find_repeats_out_merged="${find_repeats_out_root}/merged.tsv"
mock_repeats_file="${pipeline_workdir}/mock_repeats.fasta"
# TODO: remove?
# test_combinations_file="${pipeline_workdir}/test_combintations.tsv"


for dpath in "${pipeline_workdir}" "${find_repeats_out_root}"; do
    mkdir -pv "${dpath}"
done

cd "${pipeline_dir}"

# echo "$(date) -- Running extract_mock_repeats.py"
# python3 extract_mock_repeats.py \
#     "${genome_fasta}" \
#     "${mock_repeats_file}"

# echo "$(date) -- Running mutate_mock_repeats.py"
# python3 mutate_mock_repeats.py \
#     "${mock_repeats_file}" \
#     "${pipeline_workdir}"

# echo "$(date) -- Running insert_mock_repeats.py"
# python3 insert_mock_repeats.py \
#     "${genome_fasta}" \
#     "${mock_repeats_file}" \
#     "${pipeline_workdir}" \
#     "${MASK_ALL_BUT_REPEATS}"

echo "$(date) -- Running find_mock_repeats.py"
python3 find_mock_repeats_parallel.py \
    "${pipeline_workdir}" \
    "${reverlor_find_fpath}" \
    "${find_repeats_out_root}" \
    "${FINDER}" \
    "${N_THREADS}"

echo "$(date) -- Running merge_find_results.py"
python3 merge_find_results.py \
    "${find_repeats_out_root}" \
    "${find_repeats_out_merged}"

# python3 summarize_recall.py \
#     "${find_repeats_out_merged}" \
    



exit 0
