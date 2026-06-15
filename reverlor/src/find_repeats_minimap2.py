#!/usr/bin/env python3

import os
import sys
import subprocess as sp
from .FindArgs import FindArgs
from .bed_lib import merge_features


def find_repeats(args: FindArgs) -> str:

    os.makedirs(args.output_dir, exist_ok=True)

    raw_bed_fpath = os.path.join(args.output_dir, 'repeats_raw.bed')
    merged_bed_fpath = os.path.join(args.output_dir, 'repeats_final.bed')

    _self_map_genome(args, raw_bed_fpath)
    merge_features(args, raw_bed_fpath, merged_bed_fpath)

    if os.path.isfile(raw_bed_fpath):
        os.remove(raw_bed_fpath)
    # end if

    sys.stderr.write('\n')
    sys.stderr.write('INFO: Completed!\n')
    sys.stderr.write(f'INFO: output directory: `{args.output_dir}`\n')
    sys.stderr.write(f'INFO: main output file: `{merged_bed_fpath}`\n')

    return merged_bed_fpath
# end def


def _self_map_genome(args: FindArgs,
                     output_bed_fpath: str) -> None:

    cmd = [args.minimap2_fpath]
    if args.minimap_x is not None:
        cmd += ['-x', args.minimap_x]
    # end if
    cmd += [
        '-c',
        '-PD',
        '-k', str(args.minimap_k),
        '-w', str(args.minimap_w),
        '-m', str(args.minimap_m),
        '-t', '1',
        args.fasta_fpath,
        args.fasta_fpath,
    ]

    with open(output_bed_fpath, 'w') as bed_handle:
        proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
        for line in proc.stdout:
            fields = line.strip().split('\t')
            # TODO: what should I actually do in this case??
            if len(fields) < 12:
                continue
            # end if
            # Query coords (cols 1, 3, 4)
            bed_handle.write(f'{fields[0]}\t{fields[2]}\t{fields[3]}\n')
            # Target coords (cols 6, 8, 9)
            bed_handle.write(f'{fields[5]}\t{fields[7]}\t{fields[8]}\n')
        # end for
        proc.wait()
        if proc.returncode != 0:
            err = proc.stderr.read()
            sys.stderr.write(f'ERROR: minimap2 failed with code {proc.returncode}\n')
            sys.stderr.write('{}\n'.format(err))
            sys.exit(1)
        # end if
    # end with
# end def
