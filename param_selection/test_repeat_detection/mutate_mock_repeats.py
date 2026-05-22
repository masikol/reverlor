#!/usr/bin/env python3

import os
import sys

import numpy as np

from mock_repeats_settings import RATE_FROM, \
                                  RATE_TO, \
                                  RATE_STEP, \
                                  INDEL_MIN_LEN, \
                                  INDEL_MAX_LEN


infpath = os.path.abspath(sys.argv[1])
outdir_path = os.path.abspath(sys.argv[2])


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
        '--insert', rate_str,
        '--deletion', rate_str,
        '--insertminlength', str(INDEL_MIN_LEN),
        '--insertmaxlength', str(INDEL_MAX_LEN),
        '--deletionminlength', str(INDEL_MIN_LEN),
        '--deletionmaxlength', str(INDEL_MAX_LEN),
    ])

    returncode = os.system(cmd)
    if returncode != 0:
        sys.stderr.write('Error running mutation-simulator')
        sys.exit(1)
    # end if
# end def



sys.exit(0)
