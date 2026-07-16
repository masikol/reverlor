#!/usr/bin/env python3

import os
import re
import sys
import shutil
import tempfile
import subprocess as sp
from io import StringIO

from Bio import AlignIO
from Bio.Align import stockholm

from .FindArgs import FindArgs
from .bed_lib import merge_features
from .util import rm_files_if_exist


_REPEAT_ID_RE = re.compile(r'([^:]+):(\d+)-(\d+)_([+-])')
_RANDOM_SEED = 25


def find_repeats(args: FindArgs) -> str:

    os.makedirs(args.output_dir, exist_ok=True)

    raw_bed_fpath = os.path.join(args.output_dir, 'repeats_raw.bed')
    merged_bed_fpath = os.path.join(args.output_dir, 'repeats_final.bed')

    sys.stderr.write('INFO: Silently running RepeatModeler\n')
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

    remodeler_db_dir = tempfile.mkdtemp(dir=args.tmpdir)
    remodeler_result_dir = tempfile.mkdtemp(dir=args.tmpdir)

    sys.stderr.write(f'INFO: RepeatModeler database dir: `{remodeler_db_dir}`\n')
    sys.stderr.write(f'INFO: RepeatModeler work dir: `{remodeler_result_dir}`\n')

    db_name = 'blastBD'

    starting_dirpath = os.getcwd()
    os.chdir(remodeler_db_dir)
    _run_build_database(args, db_name, remodeler_db_dir)
    os.chdir(remodeler_result_dir)
    _run_repeat_modeler(args, db_name, remodeler_db_dir, remodeler_result_dir)
    os.chdir(starting_dirpath)

    rm_out_dir = _find_rm_dir(remodeler_result_dir)
    stk_fpath = os.path.join(rm_out_dir, 'families-classified.stk')

    _stk_to_bed(stk_fpath, output_bed_fpath)

    if not args.keep_tmp:
        shutil.rmtree(remodeler_db_dir, ignore_errors=True)
        shutil.rmtree(remodeler_result_dir, ignore_errors=True)
    # end if
# end def


def _run_build_database(args: FindArgs,
                        db_name: str,
                        db_dir: str) -> None:

    cmd = [
        args.build_database_fpath,
        '-name', db_name,
        args.fasta_fpath,
    ]
    proc = sp.run(cmd, cwd=db_dir, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
    if proc.returncode != 0:
        sys.stderr.write(
            'ERROR: BuildDatabase failed with code {}\n'.format(
                proc.returncode
            )
        )
        sys.stderr.write('CMD: "{}"\n'.format(' '.join(cmd)))
        sys.stderr.write(proc.stderr)
        sys.exit(1)
    # end if
# end def


def _run_repeat_modeler(args: FindArgs,
                        db_name: str,
                        db_dir: str,
                        work_dir: str) -> None:

    db_path = os.path.join(db_dir, db_name)
    cmd = [
        args.repeat_modeler_fpath,
        '-database', db_path,
        '-threads', str(args.threads),
        '-srand', str(_RANDOM_SEED),
    ]
    proc = sp.run(cmd, cwd=work_dir, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
    if proc.returncode != 0:
        sys.stderr.write(
            'ERROR: RepeatModeler failed with code {}\n'.format(
                proc.returncode
            )
        )
        sys.stderr.write('CMD: "{}"\n'.format(' '.join(cmd)))
        sys.stderr.write(proc.stderr)
        sys.exit(1)
    # end if
# end def


def _find_rm_dir(work_dir: str) -> str:

    entries = os.listdir(work_dir)
    repeatmmodeler_dirs = [
        d for d in entries
        if d.startswith('RM_') and os.path.isdir(os.path.join(work_dir, d))
    ]
    if len(repeatmmodeler_dirs) == 0:
        sys.stderr.write(
            'ERROR: no RM_* directory found in RepeatModeler work dir `{}`\n'.format(
                work_dir
            )
        )
        sys.exit(1)
    # end if
    if len(repeatmmodeler_dirs) > 1:
        sys.stderr.write(
            'ERROR: multiple ({}) RM_* directories found in RepeatModeler work dir `{}`\n'.format(
                len(repeatmmodeler_dirs),
                work_dir
            )
        )
        sys.stderr.write('Cowardly refusing to choose between them')
        sys.exit(1)
    # end if
    return os.path.join(work_dir, repeatmmodeler_dirs[0])
# end def


def _stk_to_bed(stk_fpath: str,
                output_bed_fpath: str) -> None:
    with open(output_bed_fpath, 'wt') as _:
        pass
    # end with

    if not os.path.exists(stk_fpath):
        sys.stderr.write('ERROR: final RepeatModeler file not found: `{}`\n'.format(stk_fpath))
        sys.stderr.write('Assuming that there are no repeats found')
        return
    # end if

    record_sep = '//\n'
    with open(stk_fpath, 'rt') as input_handle:
        # We will extract repeat coordinates from here
        fam_record_strings = input_handle.read().split(record_sep)
    # end with
    fam_record_strings = tuple(filter(
        lambda s: s != '',
        fam_record_strings
    ))
    if not fam_record_strings:
        sys.stderr.write('ERROR: no records in stockholm file `{}`\n'.format(stk_fpath))
        sys.exit(1)
    # end if

    # We will extract fam_id from here
    fam_records = list(stockholm.AlignmentIterator(stk_fpath))

    # TODO: proper message
    if len(fam_records) != len(fam_record_strings):
        sys.stderr.write(
            'PARSING ERROR: lengths of fam_records (len {}) and fam_record_strings (len {}) are not equal`: {}`\n' \
                .format(
                    len(fam_records),
                    len(fam_record_strings),
                    stk_fpath
                )
        )
        print('  fam_records:', file=sys.stderr)
        print(fam_records)
        print('  fam_record_strings:', file=sys.stderr)
        print(fam_record_strings)
        sys.exit(1)
    # end if

    with open(output_bed_fpath, 'wt') as bed_handle:
        for record_str, record in zip(fam_record_strings, fam_records):
            fam_id = record.annotations['identifier']
            align = AlignIO.read(StringIO(record_str), 'stockholm')
            for ar in align:
                m = _REPEAT_ID_RE.match(ar.id)
                if m is None:
                    sys.stderr.write('ERROR: cannot parse repeat_id `{}`\n'.format(ar.id))
                    sys.exit(1)
                # end if
                seqid = m.group(1)
                start = int(m.group(2)) - 1
                end   = m.group(3)
                strand = m.group(4)
                bed_handle.write('{}\n'.format('\t'.join([
                    seqid,
                    str(start),
                    end,
                    fam_id,
                    '.',
                    strand,
                ])))
            # end for
        # end for
    # end with
# end def
