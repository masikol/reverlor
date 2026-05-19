#!/usr/bin/env python3

import os
import sys
import subprocess as sp
from .FindArgs import FindArgs


def find_repeats(args: FindArgs) -> str:

    os.makedirs(args.output_dir, exist_ok=True)

    raw_bed_fpath = os.path.join(args.output_dir, 'repeats_raw.bed')
    merged_bed_fpath = os.path.join(args.output_dir, 'repeats_final.bed')

    _self_map_genome(args, raw_bed_fpath)
    _merge_repeats(raw_bed_fpath, merged_bed_fpath, args)

    # TODO: uncomment
    # if os.path.isfile(raw_bed):
    #     os.remove(raw_bed)
    # # end if

    sys.stderr.write('\n')
    sys.stderr.write('INFO: Completed!\n')
    sys.stderr.write(f'INFO: output directory: `{args.output_dir}`\n')
    sys.stderr.write(f'INFO: main output file: `{merged_bed_fpath}`\n')

    return merged_bed_fpath
# end def


def _self_map_genome(args: FindArgs,
                     output_bed_fpath: str) -> None:

    cmd = [
        args.minimap2_fpath,
        '-c',
        '-PD', '-k19', '-w19', '-m127',
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


def _merge_repeats(input_bed_fpath: str,
                   output_bed_fpath: str,
                   args: FindArgs) -> None:

    sort_cmd = [args.bedtools_fpath, 'sort', '-i', input_bed_fpath]
    merge_cmd = [args.bedtools_fpath, 'merge', '-d', '100']

    proc_sort = sp.Popen(sort_cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
    proc_merge = sp.Popen(merge_cmd, stdin=proc_sort.stdout, stdout=sp.PIPE, stderr=sp.PIPE, text=True)

    if proc_sort.stdout:
        proc_sort.stdout.close()
    # end if

    out_str, err_str = proc_merge.communicate()
    proc_sort.wait()

    if proc_sort.returncode != 0:
        sys.stderr.write('Error running bedtoools merge:\n')
        sys.stderr.write('{}\n'.format(err_str))
        sys.exit(1)
    # end if

    with open(output_bed_fpath, 'w') as out_handle:
        for line in out_str.splitlines():
            values = line.split('\t')
            start_coord = int(values[1]) # 0-based, closed
            end_coord   = int(values[2]) # 0-based, open
            region_len = end_coord - start_coord
            if region_len >= args.min_repeat_len:
                out_handle.write(line)
                out_handle.write('\n')
            # end if
        # end for
    # end with
# end def
