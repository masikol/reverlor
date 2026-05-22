#!/usr/bin/env bash

set -euo pipefail

WORKDIR='/mnt/data/Max/repos/reverlor'
reverlor_find_fpath="${WORKDIR}/reverlor/reverlor_find.py"
genome_fasta="${WORKDIR}/param_selection/data/Mycoplasma_mycoides_JCVI-syn3.0.fasta"
pipeline_dir="${WORKDIR}/param_selection/test_repeat_detection"

pipeline_data_dir="${pipeline_dir}/data"
find_repeats_out_root="${pipeline_dir}/data/find_repeats_results"
find_repeats_out_merged="${find_repeats_out_root}/merged.tsv"
mock_repeats_file="${pipeline_data_dir}/mock_repeats.fasta"
tmp_dir='/mnt/tmp_buff'


for dpath in "${pipeline_data_dir}" "${find_repeats_out_root}"; do
    mkdir -pv "${dpath}"
done

cd "${pipeline_dir}"

python3 extract_mock_repeats.py \
    "${genome_fasta}" \
    "${mock_repeats_file}"

python3 mutate_mock_repeats.py \
    "${mock_repeats_file}" \
    "${pipeline_data_dir}"

python3 insert_and_find_mock_repeats.py \
    "${genome_fasta}" \
    "${mock_repeats_file}" \
    "${pipeline_data_dir}" \
    "${reverlor_find_fpath}" \
    "${tmp_dir}" \
    "${find_repeats_out_root}"

python3 merge_find_results.py \
    "${find_repeats_out_root}" \
    "${find_repeats_out_merged}"

# python3 summarize_recall.py \
#     "${find_repeats_out_merged}" \
    



exit 0
