#!/usr/bin/env python3

import os
import sys
import glob

import numpy as np

from mock_repeats_settings import RATE_FROM, \
                                  RATE_TO, \
                                  RATE_STEP, \
                                  INDEL_MIN_LEN, \
                                  INDEL_MAX_LEN, \
                                  N_REPEAT_COPIES_TO_INSERT, \
                                  TMP_DIR_PATH


infpath = os.path.abspath(sys.argv[1])
outdir_path = os.path.abspath(sys.argv[2])
mutation_type = sys.argv[3]

assert mutation_type in ('SNP', 'SNP_INS', 'SNP_DEL', 'SNP_INDEL',)


def run_mutation_simulator(infpath: str,
                           rate_str: str,
                           mutation_type: str,
                           out_base: str) -> str:

    cmd = ' '.join([
        'mutation-simulator',
        infpath,
        '-o', out_base,
        '--quiet',
        'args',
        '--snp', rate_str,
    ])

    if mutation_type == 'SNP_INS' or mutation_type == 'SNP_INDEL':
        cmd += ' ' + ' '.join([
            '--insert', rate_str,
            '--insertminlength', str(INDEL_MIN_LEN),
            '--insertmaxlength', str(INDEL_MAX_LEN),
        ])
    # end if
    if mutation_type == 'SNP_DEL' or mutation_type == 'SNP_INDEL':
        cmd += ' ' + ' '.join([
            '--deletion', rate_str,
            '--deletionminlength', str(INDEL_MIN_LEN),
            '--deletionmaxlength', str(INDEL_MAX_LEN),
        ])
    # end if

    returncode = os.system(cmd)
    if returncode != 0:
        sys.stderr.write('Error running mutation-simulator')
        sys.exit(1)
    # end if

    rep_vcf = out_base + '_ms.vcf'
    if os.path.isfile(rep_vcf):
        os.unlink(rep_vcf)
    # end if

    out_fasta_fpath = out_base + '_ms.fasta'
    return out_fasta_fpath
# end def


for rate in np.arange(RATE_FROM, RATE_TO + RATE_STEP, RATE_STEP):
    in_base = os.path.basename(infpath).replace('.fasta', '')
    rate_str = '{:.2f}'.format(rate)
    rate_suffix = rate_str.replace('.', 'dot')

    print(rate_str)

    tmp_rep_fastas = []

    out_fasta = os.path.join(
        outdir_path,
        '{}_{}.fasta'.format(in_base, rate_suffix)
    )

    with open(out_fasta, 'wt') as out_handle:

        for copy_idx in range(1, N_REPEAT_COPIES_TO_INSERT + 1):
            rep_base = os.path.join(
                TMP_DIR_PATH,
                '{}_{}_rep{}'.format(in_base, rate_suffix, copy_idx)
            )

            tmp_rep_fasta_fpath = run_mutation_simulator(
                infpath,
                rate_str,
                mutation_type,
                rep_base
            )

            with open(tmp_rep_fasta_fpath, 'rt') as in_handle:
                for line in in_handle:
                    if line.startswith('>'):
                        line = line.replace('_c0_', '_c{}_'.format(copy_idx))
                    # end if
                    out_handle.write(line)
                # end for
                out_handle.write('\n')
            # end with

            os.unlink(tmp_rep_fasta_fpath)
        # end for
    # end with
# end for

sys.exit(0)
