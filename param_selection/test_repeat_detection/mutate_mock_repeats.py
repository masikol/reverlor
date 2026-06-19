#!/usr/bin/env python3

import os
import sys
import glob

import numpy as np

from mock_repeats_settings import RATE_FROM, \
                                  RATE_TO, \
                                  RATE_STEP, \
                                  INDEL_MIN_LEN, \
                                  INDEL_MAX_LEN


infpath = os.path.abspath(sys.argv[1])
outdir_path = os.path.abspath(sys.argv[2])
mutation_type = sys.argv[3]

assert mutation_type in ('SNP', 'SNP_INS', 'SNP_DEL', 'SNP_INDEL',)


for rate in np.arange(RATE_FROM, RATE_TO + RATE_STEP, RATE_STEP):
    in_base = os.path.basename(infpath).replace('.fasta', '')
    rate_str = '{:.2f}'.format(rate)

    print(rate_str)

    out_base = os.path.join(
        outdir_path,
        '{}_{}'.format(in_base, rate_str.replace('.', 'dot'))
    )

    print(out_base)

    cmd = ' '.join([
        'mutation-simulator',
        infpath,
        '-o', out_base,
        '--quiet',
        'args',
        '--snp', rate_str,
    ])

    if mutation_type == 'SNP_INS' or mutation_type == 'SNP_INDEL':
        cmd += ' '.join([
            '--insert', rate_str,
            '--insertminlength', str(INDEL_MIN_LEN),
            '--insertmaxlength', str(INDEL_MAX_LEN),
        ])
    # end if
    if mutation_type == 'SNP_DEL' or mutation_type == 'SNP_INDEL':
        cmd += ' '.join([
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
# end for


vcf_pattern = os.path.join(outdir_path, '*.vcf')
for fpath in glob.iglob(vcf_pattern):
    if os.path.isfile(fpath):
        os.unlink(fpath)
    # end if
# end if

sys.exit(0)
