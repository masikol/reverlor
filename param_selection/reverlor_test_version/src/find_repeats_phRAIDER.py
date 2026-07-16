#!/usr/bin/env python3

import os
import sys
import shutil
import subprocess as sp
from .FindArgs import FindArgs
from .bed_lib import merge_features
from .util import rm_files_if_exist


def find_repeats(args: FindArgs) -> str:

    os.makedirs(args.output_dir, exist_ok=True)

    raw_bed_fpath = os.path.join(args.output_dir, 'repeats_raw.bed')
    merged_bed_fpath = os.path.join(args.output_dir, 'repeats_final.bed')

    sys.stderr.write('INFO: Silently running phRAIDER\n')
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

    phraider_out_dir = os.path.join(args.output_dir, 'phraider_out')
    os.makedirs(phraider_out_dir, exist_ok=True)

    # phRAIDER with default seed is over-sensitive, as I see
    # As the authors acknowledge, it is “still not well understood what makes a good seed”
    # https://github.com/karroje/phRAIDER
    seed = '1' * args.min_repeat_len

    phraider_cmd = [
        args.phraider_fpath,
        '-s', seed,
        '-m', str(args.min_repeat_len),
        '-c', '2', # find two repeats per family or more
        args.fasta_fpath,
        phraider_out_dir,
    ]

    sys.stderr.write('Runnging command:\n  `{}`\n'.format(' '.join(phraider_cmd)))

    phraider_proc = sp.run(
        phraider_cmd,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        text=True
    )
    if phraider_proc.returncode != 0:
        sys.stderr.write(
            'ERROR: phRAIDER failed with code {}\n'.format(
                phraider_proc.returncode
            )
        )
        sys.stderr.write('CMD: “{}”\n'.format(' '.join(phraider_cmd)))
        sys.stderr.write(phraider_proc.stderr)
        sys.exit(1)
    # end if

    # Convert elements file to BED
    elements_fpath = os.path.join(phraider_out_dir, 'elements')
    _elements_to_bed(elements_fpath, output_bed_fpath)

    if not args.keep_tmp:
        shutil.rmtree(phraider_out_dir)
    # end if
# end def


def _elements_to_bed(elements_fpath: str,
                     output_bed_fpath: str) -> None:
    sep = '\t'

    # elements columns:
    #   1. repeat_family -> BED col 4 (name)
    #   2. element -> ignore
    #   3. direction -> always "1", ignore
    #   4. sequence name -> first word  -> BED col 1 (chrom)
    #   5. start -> BED col 2 (0-based, closed, unsure)
    #   6. end   -> BED col 3 (0-based, open, unsure)
    with open(elements_fpath, 'rt') as in_handle, \
         open(output_bed_fpath, 'wt') as bed_handle:

        for line in in_handle:
            if line.startswith('#'):
                continue
            # end if

            fields = line.strip().split(sep)

            if len(fields) < 6:
                sys.stderr.write(
                    'ERROR: corrupted line in the elements file `{}`\n'.format(
                        elements_fpath
                    )
                )
                sys.stderr.write('The line: “{}”\n'.format(line))
                sys.exit(1)
            # end if

            family     = fields[0]
            seq_full   = fields[3]
            chrom      = seq_full.partition(' ')[0]
            start      = fields[4]
            end        = fields[5]

            bed_handle.write('{}\n'.format(sep.join([
                chrom,
                start,
                end,
                'fam:{}'.format(family),
                '.',
                '.',
            ])))
        # end for
    # end with
# end def
