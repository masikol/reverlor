#!/usr/bin/env python3

import os
import re
import sys
import glob

import polars as pl


in_detected_dir = os.path.abspath(sys.argv[1])
in_true_dir = os.path.abspath(sys.argv[2])
out_fpath = os.path.abspath(sys.argv[3])


IN_GLOB = sorted(glob.glob(
    os.path.join(in_true_dir, 'MR_l*_r*_t*.bed')
))

FIND_RESULT_DIR_NAME_RE = re.compile(
    r'^MR_l(\d+)_r(\d+)_t([0-9\.]+)$'
)


sep = '\t'

with open(out_fpath, 'wt') as out_handle:

    out_handle.write(
        '{}\n'.format(sep.join([
            'chr',
            'repeat_len',
            'repeat_idx',
            'rate',
            'value_type',
            'start_coord',
            'end_coord',
            'strand',
        ]))
    )

    for true_fpath in IN_GLOB:
        replicate_id = os.path.basename(true_fpath.replace('.bed', ''))
        re_obj = re.match(
            FIND_RESULT_DIR_NAME_RE,
            replicate_id
        )

        repeat_len = str(re_obj.group(1))
        replicate_idx = str(re_obj.group(2))
        rate = str(re_obj.group(3))

        # TODO: remove
        # dir_base_name = os.path.basename(in_dir_path)
        # true_fpath = os.path.join(in_true_dir, dir_base_name, 'repeats_final.true.bed')
        pred_fpath = os.path.join(in_detected_dir, replicate_id, 'repeats_final.bed')
        out_zip = zip(
            [true_fpath, pred_fpath,],
            ['true', 'pred',],
        )

        for bed_fpath, value_type in out_zip:

            with open(bed_fpath, 'rt') as in_handle:
                lines = in_handle.readlines()
            # end def

            for line in lines:
                bed_values = line.strip().split('\t')

                strand = '.'
                if len(bed_values) > 4:
                    strand = bed_values[4]
                # end if

                out_handle.write('{}\n'.format(sep.join([
                    bed_values[0],
                    repeat_len,
                    replicate_idx,
                    rate,
                    value_type,
                    bed_values[1],
                    bed_values[2],
                    strand,
                ])))
            # end for
        # end for

    # end for
# end with


print('Completed!')
print(out_fpath)
sys.exit(0)
