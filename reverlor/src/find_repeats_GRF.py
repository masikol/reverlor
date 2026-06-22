#!/usr/bin/env python3

import os
import sys
import shutil
import tempfile
import subprocess as sp

from .FindArgs import FindArgs
from .bed_lib import merge_features
from .util import rm_files_if_exist


def find_repeats(args: FindArgs) -> str:

    os.makedirs(args.output_dir, exist_ok=True)

    raw_bed_fpath = os.path.join(args.output_dir, 'repeats_raw.bed')
    merged_bed_fpath = os.path.join(args.output_dir, 'repeats_final.bed')

    sys.stderr.write('INFO: Silently running GRF\n')
    _create_raw_repeat_file(args, raw_bed_fpath)
    sys.stderr.write('INFO: Silently merging repeats\n')
    merge_features(args, raw_bed_fpath, merged_bed_fpath)

    if not args.keep_tmp:
        rm_files_if_exist(raw_bed_fpath)
    # end if

    sys.stderr.write('\n')
    sys.stderr.write('INFO: Completed!\n')
    sys.stderr.write(f'INFO: output directory: `{args.output_dir}`\n')
    sys.stderr.write(f'INFO: main output file: `{merged_bed_fpath}`\n')

    return merged_bed_fpath
# end def


def _create_raw_repeat_file(args: FindArgs,
                            output_bed_fpath: str) -> None:

    tmp_dir = tempfile.mkdtemp(dir=args.tmpdir)

    out_fmt_mode = '1'
    min_seed_match_num = '2'
    seed_region_len = '100' # bp

    grf_cmd = [
        args.grf_intersperse_fpath,
        '-i', args.fasta_fpath,
        '-o', tmp_dir,
        '-f', out_fmt_mode,
        '-c', min_seed_match_num,
        '-t', str(args.threads),
        '-s', seed_region_len,
    ]
    grf_proc = sp.run(grf_cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
    if grf_proc.returncode != 0:
        sys.stderr.write(
            'ERROR: grf-intersperse failed with code {}\n'.format(
                grf_proc.returncode
            )
        )
        sys.stderr.write('CMD: “{}”\n'.format(' '.join(grf_cmd)))
        sys.stderr.write(grf_proc.stderr)
        sys.exit(1)
    # end if

    out_fpath = os.path.join(tmp_dir, 'interspersed_repeat.out')
    _out_to_bed(out_fpath, output_bed_fpath)

    if not args.keep_tmp:
        shutil.rmtree(tmp_dir, ignore_errors=True)
    # end if
# end def


def _out_to_bed(grf_out_fpath: str,
                output_bed_fpath: str) -> None:
    sep = '\t'

    # .out lines starting with '>'
    # Format: >chrom:start:end:strand
    # Example: >CP002027.1:1002259:1002361:-
    with open(grf_out_fpath, 'rt') as in_handle, \
         open(output_bed_fpath, 'wt') as bed_handle:

        for line in in_handle:
            if not line.startswith('>'):
                continue
            # end if

            fields = line[1:].strip().split(':')
            if len(fields) < 4:
                sys.stderr.write(
                    'PARSING ERROR: cannot parse line in file `{}`\n'.format(
                        grf_out_fpath
                    )
                )
                sys.stderr.write('Here is the line: `{}`\n'.format(line))
                sys.exit(1)
            # end if

            chrom  = fields[0]
            start  = int(fields[1]) - 1  # 1-based closed -> 0-based closed
            end    = int(fields[2])       # 1-based closed -> 0-based open
            strand = fields[3]

            bed_handle.write('{}\n'.format(sep.join([
                chrom,
                str(start),
                str(end),
                '.',
                '.',
                strand,
            ])))
        # end for
    # end with
# end def
