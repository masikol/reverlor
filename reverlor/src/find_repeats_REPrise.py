#!/usr/bin/env python3

import os
import sys
import tempfile
import subprocess as sp

from .FindArgs import FindArgs
from .bed_lib import merge_features
from .util import rm_files_if_exist, rm_empty_dir_if_exists


def find_repeats(args: FindArgs) -> str:

    os.makedirs(args.output_dir, exist_ok=True)
    tmp_dir_path = tempfile.mkdtemp(dir=args.tmpdir)

    sys.stderr.write('INFO: Silently running REPrise\n')
    _create_raw_repeat_file(args, tmp_dir_path)
    sys.stderr.write('INFO: Silently merging repeats\n')

    raw_bed_fpath = os.path.join(tmp_dir_path, 'REPrise_result.bed')
    _fix_first_bed_column(args, raw_bed_fpath)
    merged_bed_fpath = os.path.join(args.output_dir, 'repeats_final.bed')
    merge_features(args, raw_bed_fpath, merged_bed_fpath)

    if not args.keep_tmp:
        rm_files_if_exist(raw_bed_fpath)
        rm_empty_dir_if_exists(tmp_dir_path)
    # end if

    sys.stderr.write('\n')
    sys.stderr.write('INFO: Completed!\n')
    sys.stderr.write(f'INFO: output directory: `{args.output_dir}`\n')
    sys.stderr.write(f'INFO: main output file: `{merged_bed_fpath}`\n')

    return merged_bed_fpath
# end def


def _create_raw_repeat_file(args: FindArgs,
                            outdir: str) -> None:

    out_prefix = os.path.join(outdir, 'REPrise_result')

    reprise_cmd = [
        args.reprise_fpath,
        '-input', args.fasta_fpath,
        '-output', out_prefix,
        '-additonalfile',
        '-minlength', str(args.min_repeat_len),
        '-minfreq', '2',
        '-pa', str(args.threads),
    ]

    sys.stderr.write('Runnging command:\n  `{}`\n'.format(' '.join(reprise_cmd)))

    reprise_proc = sp.run(
        reprise_cmd,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        text=True,
    )
    if reprise_proc.returncode != 0:
        sys.stderr.write(
            'ERROR: REPrise failed with code {}\n'.format(
                reprise_proc.returncode
            )
        )
        sys.stderr.write('CMD: “{}”\n'.format(' '.join(reprise_cmd)))
        sys.stderr.write(reprise_proc.stderr)
        sys.exit(1)
    # end if

    if not args.keep_tmp:
        for ext in ('.freq', '.masked', '.reprof'):
            rm_files_if_exist('{}{}'.format(out_prefix, ext))
        # end for
    # end if
# end def


def _fix_first_bed_column(args: FindArgs,
                          raw_bed_fpath: str) -> None:

    # collect sequence IDs from the FASTA file
    seq_ids = []
    with open(args.fasta_fpath, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                seq_id = line[1:].strip().split()[0]
                seq_ids.append(seq_id)
            # end if
        # end for
    # end with

    seq_ids.sort(key=len, reverse=True)  # longer first to avoid prefix collisions

    raw_bed_tmp = raw_bed_fpath + '.tmp'
    with open(raw_bed_fpath, 'rt') as in_handle, \
         open(raw_bed_tmp, 'wt') as out_handle:

        for line in in_handle:
            fields = line.rstrip('\n').split('\t')

            seq_id = _find_matching_seq_id(fields[0], seq_ids)
            fields[0] = seq_id

            out_handle.write('\t'.join(fields) + '\n')
        # end for
    # end with

    os.replace(raw_bed_tmp, raw_bed_fpath)
# end def


def _find_matching_seq_id(chrom_value: str, seq_ids: list[str]) -> str:
    for seq_id in seq_ids:
        if chrom_value.startswith(f'{seq_id}_'):
            return seq_id
        # end if
    # end for
    sys.stderr.write(
        'ERROR: cannot match raw BED col1 "{}" to any input FASTA sequence ID\n'
        .format(chrom_value)
    )
    sys.stderr.write(
        'Here are the detected input FASTA sequence IDs: {}\n'
        .format(str(seq_ids))
    )
    sys.exit(1)
# end def
