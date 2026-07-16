#!/usr/bin/env python3

import os
import re
import sys
import glob
import tempfile
import subprocess as sp

from .FindArgs import FindArgs
from .bed_lib import merge_features
from .util import rm_files_if_exist, rm_empty_dir_if_exists


def find_repeats(args: FindArgs) -> str:

    os.makedirs(args.output_dir, exist_ok=True)

    raw_bed_fpath = os.path.join(args.output_dir, 'repeats_raw.bed')
    merged_bed_fpath = os.path.join(args.output_dir, 'repeats_final.bed')

    sys.stderr.write('INFO: Silently running TotalRepeats\n')
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

    tmp_dir_extract = tempfile.mkdtemp(dir=args.tmpdir)
    tmp_dir_repeats = tempfile.mkdtemp(dir=args.tmpdir)

    _run_extract_fasta(args, tmp_dir_extract)
    _run_detect_repeats(args, tmp_dir_extract, tmp_dir_repeats)

    # Convert GFF files to BED
    gff_fpaths = _find_relevant_gff_files(tmp_dir_repeats)

    _gffs_to_bed(gff_fpaths, output_bed_fpath)

    if not args.keep_tmp:
        _clean_up(tmp_dir_extract, tmp_dir_repeats)
    # end if
# end def


def _run_extract_fasta(args: FindArgs,
                       tmp_dir: str) -> None:

    extract_cmd = [
        args.java_fpath, '-jar', args.total_repeats_fpath,
        args.fasta_fpath,
        '-out={}'.format(tmp_dir),
        '-extract',
    ]

    if args.threads == 1:
        extract_cmd.append('normal')
    # end if

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
# end def


def _run_detect_repeats(args: FindArgs,
                        tmp_dir_extract: str,
                        tmp_dir_repeats: str) -> None:

    detect_cmd = [
        args.java_fpath, '-jar', args.total_repeats_fpath,
        tmp_dir_extract,
        '-out={}'.format(tmp_dir_repeats),
        '-sln={}'.format(args.min_repeat_len),
        '-joint',
    ]

    if args.threads == 1:
        detect_cmd.append('normal')
    # end if

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
# end def


def _find_relevant_gff_files(dir_path: str) -> list[str]:
    gff_pattern = os.path.join(dir_path, '*.fasta.gff')
    gff_fpaths = sorted(glob.glob(gff_pattern))
    if len(gff_fpaths) == 0:
        sys.stderr.write(
            'ERROR: no *.fasta.gff files found in `{}`\n'.format(
                args.output_dir
            )
        )
        sys.exit(1)
    # end if
    return gff_fpaths
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
                # Pass first 12 lines
                for _ in range(12):
                    gff_handle.readline()
                # end for

                for line_idx, line in enumerate(gff_handle, 1):

                    fields = line.strip().split(sep)

                    if len(fields) < 7:
                        sys.stderr.write(
                            'ERROR: cannot parse line #{} in file `{}`\n'.format(
                                line_idx,
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


def _clean_up(tmp_dir_extract: str, tmp_dir_repeats: str) -> None:
    fpaths_to_rm = glob.glob(os.path.join(tmp_dir_extract, '*.fasta')) \
        + glob.glob(os.path.join(tmp_dir_extract, 'report.txt'))
    for ext in ('*.gff', '*.svg', '*.png', '*.msk', '*.tsv', '*.txt'):
        fpaths_to_rm.extend(glob.glob(os.path.join(tmp_dir_repeats, ext)))
    # end for
    rm_files_if_exist(*fpaths_to_rm)
    rm_empty_dir_if_exists(tmp_dir_extract)
    rm_empty_dir_if_exists(tmp_dir_repeats)
# end def
