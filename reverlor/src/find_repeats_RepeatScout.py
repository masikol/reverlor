#!/usr/bin/env python3

import os
import sys
import subprocess as sp
from .FindArgs import FindArgs
from .bed_lib import merge_features
from .util import rm_files_if_exist


def find_repeats(args: FindArgs) -> str:

    os.makedirs(args.output_dir, exist_ok=True)

    raw_bed_fpath = os.path.join(args.output_dir, 'repeats_raw.bed')
    merged_bed_fpath = os.path.join(args.output_dir, 'repeats_final.bed')

    sys.stderr.write('INFO: Silently running RepearScout\n')
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

    basename = os.path.splitext(os.path.basename(args.fasta_fpath))[0]
    freq_fpath = os.path.join(args.output_dir, f'{basename}.freq')
    fasta_out_fpath = os.path.join(args.output_dir, f'{basename}.out')
    range_fpath = os.path.join(args.output_dir, f'{basename}.range')

    # Step 1: build_lmer_table
    _run_build_lmer_table(args, freq_fpath)

    # Step 2: RepeatScout
    _run_repeat_scout(args, freq_fpath, fasta_out_fpath, range_fpath)

    if not args.keep_tmp:
        rm_files_if_exist(freq_fpath, fasta_out_fpath)
    # end if

    # Step 3: Convert .range file to a BED file
    _range_to_bed(range_fpath, output_bed_fpath)

    if not args.keep_tmp:
        rm_files_if_exist(range_fpath)
    # end if
# end def

def _run_build_lmer_table(args: FindArgs, freq_fpath: str) -> None:
    build_lmer_table_cmd = [
        args.build_lmer_table_fpath,
        '-sequence', args.fasta_fpath,
        '-freq', freq_fpath,
        '-min', '2',
    ]
    build_lmer_table_proc = sp.run(
        build_lmer_table_cmd,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        text=True
    )
    if build_lmer_table_proc.returncode != 0:
        sys.stderr.write(
            'ERROR: build_lmer_table failed with code {}\n'.format(
                build_lmer_table_proc.returncode
            )
        )
        sys.stderr.write('CMD: “{}”\n'.format(' '.join(build_lmer_table_cmd)))
        sys.stderr.write(build_lmer_table_proc.stderr)
        sys.exit(1)
    # end if
# end def

def _run_repeat_scout(args: FindArgs,
                      freq_fpath: str,
                      fasta_out_fpath: str,
                      range_fpath: str) -> None:
    repeat_scout_cmd = [
        args.repeat_scout_fpath,
        '-sequence', args.fasta_fpath,
        '-freq', freq_fpath,
        '-output', fasta_out_fpath,
        '-ranges', range_fpath,
        '-goodlength', str(args.min_repeat_len),
    ]
    repeat_scout_proc = sp.run(
        repeat_scout_cmd,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        text=True
    )
    if repeat_scout_proc.returncode != 0:
        sys.stderr.write(
            'ERROR: RepeatScout failed with code {}\n'.format(
                repeat_scout_proc.returncode
            )
        )
        sys.stderr.write('CMD: {}\n'.format(' '.join(repeat_scout_cmd)))
        sys.stderr.write(repeat_scout_proc.stderr)
        sys.exit(1)
    # end if
# end def

def _range_to_bed(range_fpath: str,
                  output_bed_fpath: str) -> None:
    sep = '\t'

    # .range columns:
    #   1. Sequence ID
    #   2. Start (1-based, closed)
    #   3. End   (1-based, closed)
    #   4. Repeat family number
    #   5. Strand
    #   6-8. (don't know what they are, ignored)
    # BED conversion: start-1, end, family -> name, strand -> strand
    with open(range_fpath, 'rt') as range_handle, \
         open(output_bed_fpath, 'wt') as bed_handle:

        for line in range_handle:
            fields = line.strip().split(sep)

            if len(fields) < 5:
                sys.stderr.write(
                    'Error: invalid line in the .range file `{}`\n'.format(range_fpath)
                )
                sys.stderr.write(line)
                sys.exit(1)
            # end if

            seq_id = fields[0]
            start  = int(fields[1]) - 1  # 1-based closed -> 0-based
            end    = int(fields[2])      # closed -> open (half-open)
            family = fields[3]
            strand = fields[4]
            bed_handle.write('{}\n'.format(sep.join([
                seq_id,
                str(start),
                str(end),
                family,
                '.',
                strand
            ])))
        # end for
    # end with
# end def
