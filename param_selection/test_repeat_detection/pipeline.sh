#!/usr/bin/env bash

set -euo pipefail


# >>> Configuration >>>

FINDERS=(
    minimap2
    # repeat-scout
    # phraider
    # total-repeats
    # repeat-modeler
    # grf
    # 'reprise'
)

MUTATION_TYPE='SNP' # choice
# MUTATION_TYPE='SNP_INS' # choice
# MUTATION_TYPE='SNP_DEL' # choice
# MUTATION_TYPE='SNP_INDEL' # choice

N_DESIRED_REPEAT_COPIES=3

N_THREADS=6

WORKDIR='/mnt/data/Max/repos/reverlor'
reverlor_find_fpath="${WORKDIR}/reverlor/reverlor_find.py"

# <<< Configuration <<<


genome_fasta="${WORKDIR}/param_selection/data/Mycoplasma_mycoides_JCVI-syn3.0.fasta"
plasmid_fasta="${WORKDIR}/param_selection/data/pUC18.fasta"

pipeline_dir="${WORKDIR}/param_selection/test_repeat_detection"
mock_repeats_dir="${pipeline_dir}/mock_repeats_${MUTATION_TYPE}_${N_DESIRED_REPEAT_COPIES}copies"
workdir_root="${pipeline_dir}/workdirs_${N_DESIRED_REPEAT_COPIES}copies"
mock_repeats_file="${mock_repeats_dir}/mock_repeats.fasta"
true_repeat_dir="${mock_repeats_dir}/true_repeat_locations"
replicate_id_list="${mock_repeats_dir}/replicate_id_list.txt"

mkdir -pv "${mock_repeats_dir}"


cd "${pipeline_dir}"

echo "$(date) -- Running extract_mock_repeats.py"
python3 extract_mock_repeats.py \
    "${genome_fasta}" \
    "${mock_repeats_file}" \
    "${N_DESIRED_REPEAT_COPIES}"

echo "$(date) -- Running mutate_mock_repeats.py"
python3 mutate_mock_repeats.py \
    "${mock_repeats_file}" \
    "${mock_repeats_dir}" \
    "${MUTATION_TYPE}"

echo "$(date) -- Running insert_mock_repeats.py"
python3 insert_mock_repeats.py \
    "${genome_fasta}" \
    "${plasmid_fasta}" \
    "${mock_repeats_dir}" \
    "${N_DESIRED_REPEAT_COPIES}"


for finder in "${FINDERS[@]}"; do

    pipeline_workdir="${workdir_root}/workdir_${finder}_${MUTATION_TYPE}"
    detected_repeat_dir="${pipeline_workdir}/detected_repeat_locations"
    find_repeats_out_merged="${pipeline_workdir}/repeat_detection_table_raw.tsv"
    repeat_detection_table="${pipeline_workdir}/repeat_detection_table.tsv"

    for dpath in "${pipeline_workdir}" "${detected_repeat_dir}"; do
        mkdir -pv "${dpath}"
    done

    echo "  >>> FINDER: ${finder} >>>"

    echo "$(date) -- Running find_mock_repeats.py"
    python3 find_mock_repeats.py \
        "${mock_repeats_dir}" \
        "${pipeline_workdir}" \
        "${replicate_id_list}" \
        "${reverlor_find_fpath}" \
        "${detected_repeat_dir}" \
        "${finder}" \
        "${N_THREADS}"

    echo "$(date) -- Running merge_find_results.py"
    python3 merge_find_results.py \
        "${detected_repeat_dir}" \
        "${true_repeat_dir}" \
        "${find_repeats_out_merged}"

    echo "$(date) -- Running make_repeat_detection_table.py"
    python3 make_repeat_detection_table.py \
        "${find_repeats_out_merged}" \
        "${repeat_detection_table}" \
        "${finder}"

    echo "  <<< FINDER: ${finder} <<<"
done



exit 0
