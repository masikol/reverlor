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
