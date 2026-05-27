#!/usr/bin/env python3

import os
import re
import sys
import glob

import polars as pl


input_dir = os.path.abspath(sys.argv[1])
out_fpath = os.path.abspath(sys.argv[2])


IN_GLOB = sorted(glob.glob(
    os.path.join(input_dir, 'default_mock_repeat_*')
))

FIND_RESULT_DIR_NAME_RE = re.compile(
    r'^default_mock_repeat_([1-9][0-9]*)_([0-9]+)[_]?(0\.[0-9]+)?$'
)


sep = '\t'

with open(out_fpath, 'wt') as out_handle:

    out_handle.write(
        '{}\n'.format(sep.join([
            'repeat_len',
            'repeat_idx',
            'rate',
            'value_type',
            'start_coord',
            'end_coord',
        ]))
    )

    for in_dir_path in IN_GLOB:
        re_obj = re.match(
            FIND_RESULT_DIR_NAME_RE,
            os.path.basename(in_dir_path)
        )

        repeat_len = str(re_obj.group(1))
        repeat_idx = str(re_obj.group(2))
        rate = str(re_obj.group(3))
        if rate == 'None':
            rate = '0.00'
        # end if

        true_fpath = os.path.join(in_dir_path, 'repeats_final.true.bed')
        pred_fpath = os.path.join(in_dir_path, 'repeats_final.bed')
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
                out_handle.write('{}\n'.format(sep.join([
                    repeat_len,
                    repeat_idx,
                    rate,
                    value_type,
                    bed_values[1],
                    bed_values[2],
                ])))
            # end for
        # end for

    # end for
# end with


print('Completed!')
print(out_fpath)
sys.exit(0)
