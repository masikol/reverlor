#!/usr/bin/env python3

import os
import re
import sys
import glob
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

    sys.stderr.write('INFO: Silently running TotalRepeats\n')
    _find_repeats_total_repeats(args, raw_bed_fpath)
    sys.stderr.write('INFO: Silently merging repeats\n')
    merge_features(args, raw_bed_fpath, merged_bed_fpath)

    rm_files_if_exist(raw_bed_fpath)

    sys.stderr.write('\n')
    sys.stderr.write('INFO: Completed!\n')
    sys.stderr.write(f'INFO: output directory: `{args.output_dir}`\n')
    sys.stderr.write(f'INFO: main output file: `{merged_bed_fpath}`\n')

    return merged_bed_fpath
# end def


def _find_repeats_total_repeats(args: FindArgs,
                                output_bed_fpath: str) -> None:

    tmp_dir = tempfile.mkdtemp()

    # Step 1: extract
    extract_cmd = [
        args.java_fpath, '-jar', args.total_repeats_fpath,
        args.fasta_fpath,
        '-out={}'.format(tmp_dir),
        '-extract',
    ]
    extract_proc = sp.run(extract_cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
    if extract_proc.returncode != 0:
        sys.stderr.write(
            'ERROR: TotalRepeats (extract) failed with code {}\n'.format(
                extract_proc.returncode
            )
        )
        sys.stderr.write('CMD: “{}”\n'.format(' '.join(extract_cmd)))
        sys.stderr.write(extract_proc.stderr)
        sys.exit(1)
    # end if

    # Step 2: detect repeats
    detect_cmd = [
        args.java_fpath, '-jar', args.total_repeats_fpath,
        tmp_dir,
        '-out={}'.format(args.output_dir),
        '-sln={}'.format(args.min_repeat_len),
        '-joint',
    ]
    detect_proc = sp.run(detect_cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
    if detect_proc.returncode != 0:
        sys.stderr.write(
            'ERROR: TotalRepeats (-joint) failed with code {}\n'.format(
                detect_proc.returncode
            )
        )
        sys.stderr.write('CMD: “{}”\n'.format(' '.join(detect_cmd)))
        sys.stderr.write(detect_proc.stderr)
        sys.exit(1)
    # end if

    # Convert GFF files to BED
    gff_pattern = os.path.join(args.output_dir, '*.fasta.gff')
    gff_fpaths = sorted(glob.glob(gff_pattern))
    if len(gff_fpaths) == 0:
        sys.stderr.write(
            'ERROR: no *.fasta.gff files found in `{}`\n'.format(
                args.output_dir
            )
        )
        sys.exit(1)
    # end if
    _gffs_to_bed(gff_fpaths, output_bed_fpath)

    # Clean up files produced by TotalRepeats
    for ext in ('*.gff', '*.svg', '*.png', '*.msk', '*.tsv', '*.txt'):
        for fpath in glob.glob(os.path.join(args.output_dir, ext)):
            rm_files_if_exist(fpath)
        # end for
    # end for

    shutil.rmtree(tmp_dir, ignore_errors=True)
# end def


def _gffs_to_bed(gff_fpaths: list,
                 output_bed_fpath: str) -> None:
    sep = '\t'
    allowed_repeat_types = frozenset(['CRP', 'UCRP'])

    # GFF columns:
    #   1. sequence name -> BED col 1 (chrom)
    #   2. repeat type   -> only 'CRP' or 'UCRP' (skip otherwise)
    #   3. family ID     -> BED col 4 as "fam:{id}"
    #   4. start         -> 1-based closed -> BED col 2 (0-based, closed): start - 1
    #   5. end           -> 1-based closed -> BED col 3 (0-based, open)  : end
    #   6. ignore
    #   7. strand        -> BED col 6
    #   8. ignore
    with open(output_bed_fpath, 'wt') as bed_handle:
        for gff_fpath in gff_fpaths:
            with open(gff_fpath, 'rt') as gff_handle:
                for idx, line in enumerate(gff_handle):
                    if idx < 12:
                        continue
                    # end if

                    fields = line.strip().split(sep)

                    if len(fields) < 7:
                        sys.stderr.write(
                            'ERROR: cannot parse line #{} in file `{}`\n'.format(
                                idx + 1,
                                gff_fpath
                            )
                        )
                        sys.stderr.write(
                            'Expected >=7 tab-sparated fields, found {}\n'.format(len(fields))
                        )
                        sys.exit(1)
                    # end if

                    repeat_type = fields[1]
                    if repeat_type not in allowed_repeat_types:
                        continue
                    # end if

                    chrom       = fields[0]
                    family_id   = fields[2]
                    start       = int(fields[3]) - 1  # 1-based closed -> 0-based closed
                    end         = int(fields[4])       # 1-based closed -> 0-based open
                    strand      = fields[6]

                    bed_handle.write('{}\n'.format(sep.join([
                        chrom,
                        str(start),
                        str(end),
                        'fam:{}'.format(family_id),
                        '.',
                        strand,
                    ])))
                # end for
            # end with
        # end for
    # end with
# end def
