#!/usr/bin/env python3

import os
import sys
import subprocess as sp

import polars as pl

from mock_repeats_settings import REPEAT_SCOUT_DIR_PATH, PHRAIDER_PATH


PIPELINE_DATA_DIRPATH = os.path.abspath(sys.argv[1])
REVERLOR_FIND_FPATH = os.path.abspath(sys.argv[2])
FIND_REPEATS_OUT_ROOT_PATH = os.path.abspath(sys.argv[3])
FINDER = sys.argv[4]


def run_reverlor_find(input_fasta_fpath, finder, out_dir_path):
    cmd = [
        'python3', REVERLOR_FIND_FPATH,
        '--finder', FINDER,
        '--minimap-m', '65',
        '--min-repeat-len', '127',
        input_fasta_fpath,
        out_dir_path,
    ]

    if FINDER == 'phraider':
        cmd += ['--phraider', PHRAIDER_PATH]
    elif FINDER == 'repeat-scout':
        cmd += ['--repeat-scout', REPEAT_SCOUT_DIR_PATH]
    # end if

    pipe = sp.Popen(cmd, text=True, stdout=sp.PIPE, stderr=sp.PIPE)

    _, stderr = pipe.communicate()
    if pipe.returncode != 0:
        sys.stderr.write('Error running reverlor_find.py')
        sys.stderr.write(stderr)
        sys.exit(1)
    # end if

    return os.path.join(out_dir_path, 'repeats_final.bed')
# end def

def make_curr_outdir_path(repeat_record_id, rate=None):
    if rate == 0.0:
        base_name = 'default_{}'.format(repeat_record_id)
    else:
        base_name = 'default_{}_{:.2f}'.format(repeat_record_id, rate)
    # end def
    return os.path.join(
        FIND_REPEATS_OUT_ROOT_PATH,
        base_name
    )
# end def



test_combinations_fpath = os.path.join(
    PIPELINE_DATA_DIRPATH,
    'test_combintations.tsv'
)
chr_with_inserts_dir_path= os.path.join(
    PIPELINE_DATA_DIRPATH,
    'chr_with_inserts'
)
comb_df = pl.read_csv(test_combinations_fpath, separator='\t')

for row in comb_df.to_dicts():

    input_fasta_fpath = os.path.join(
        chr_with_inserts_dir_path,
        row['chr_with_insert_fname']
    )
    out_dir_path = make_curr_outdir_path(
        row['repeat_id'],
        row['rate']
    )

    print(row['repeat_id'], row['rate'])

    run_reverlor_find(input_fasta_fpath, FINDER, out_dir_path)
# end for


sys.exit(0)
