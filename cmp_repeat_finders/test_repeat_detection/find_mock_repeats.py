#!/usr/bin/env python3

import os
import sys
import subprocess as sp
import multiprocessing

import polars as pl

from mock_repeats_settings import MINIMAP2_FPATH, \
                                  REPEAT_SCOUT_DIR_PATH, \
                                  PHRAIDER_PATH, \
                                  TORALREPEATS_PATH, \
                                  REPEATMODELER_DIR_PATH, \
                                  BEDTOOLS_PATH, \
                                  GRF_INTERSPERSE_FPATH, \
                                  REPRISE_FPATH, \
                                  TMP_DIR_PATH


IN_MOCK_REPEATS_DIRPATH = os.path.abspath(sys.argv[1])
PIPELINE_DATA_DIRPATH = os.path.abspath(sys.argv[2])
REPLICATE_ID_LIST_FPATH = os.path.abspath(sys.argv[3])
REVERLOR_FIND_FPATH = os.path.abspath(sys.argv[4])
FIND_REPEATS_OUT_ROOT_PATH = os.path.abspath(sys.argv[5])
FINDER = sys.argv[6]
NUM_THREADS = int(sys.argv[7]) if len(sys.argv) > 6 else 1


def run_reverlor_find(input_fasta_fpath, finder, out_dir_path):
    cmd = [
        'python3', REVERLOR_FIND_FPATH,
        '--finder', FINDER,
        '--min-repeat-len', '127',
        '--tmpdir', TMP_DIR_PATH,
        '--bedtools', BEDTOOLS_PATH,
        '--threads', '1',
        input_fasta_fpath,
        out_dir_path,
    ]

    if FINDER == 'minimap2':
        cmd += ['--minimap2', MINIMAP2_FPATH,]
        cmd += ['--minimap-m', '65',]
    elif FINDER == 'phraider':
        cmd += ['--phraider', PHRAIDER_PATH]
    elif FINDER == 'repeat-scout':
        cmd += ['--repeat-scout', REPEAT_SCOUT_DIR_PATH]
    elif FINDER == 'total-repeats':
        cmd += ['--total-repeats', TORALREPEATS_PATH]
    elif FINDER == 'repeat-modeler':
        cmd += ['--repeat-modeler', REPEATMODELER_DIR_PATH]
    elif FINDER == 'grf':
        cmd += ['--grf-intersperse', GRF_INTERSPERSE_FPATH]
    elif FINDER == 'reprise':
        cmd += ['--reprise', REPRISE_FPATH]
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

def make_curr_outdir_path(replicate_id):
    return os.path.join(
        FIND_REPEATS_OUT_ROOT_PATH,
        replicate_id
    )
# end def


def process_replicate(replicate_id):
    try:
        input_fasta_fpath = os.path.join(
            chr_with_inserts_dir_path,
            '{}.fasta'.format(replicate_id)
        )
        out_dir_path = make_curr_outdir_path(
            replicate_id
        )

        print(replicate_id)

        return run_reverlor_find(input_fasta_fpath, FINDER, out_dir_path)
    except Exception as exc:
        return f'ERROR: {replicate_id}: {exc}'
    # end try
# end def



chr_with_inserts_dir_path = os.path.join(
    IN_MOCK_REPEATS_DIRPATH,
    'chr_with_inserts'
)
with open(REPLICATE_ID_LIST_FPATH, 'rt') as in_handle:
    all_replicate_ids = [
        line.strip() for line in in_handle
    ]
# end with

with multiprocessing.Pool(processes=NUM_THREADS) as pool:
    results = pool.map(process_replicate, all_replicate_ids)
# end with

# Reap zombies
pool.close()
pool.join()

for res in results:
    if isinstance(res, str) and res.startswith('ERROR:'):
        sys.stderr.write(res + '\n')
        sys.exit(1)
    # end if
# end for


sys.exit(0)
