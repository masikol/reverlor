#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess as sp
from typing import Optional


SAMTOOLS_DEFAULT_FPATH = 'samtools'
DEFAULT_NUM_READ_THRESHOLD = 5
DEFAULT_SHOULDER_LEN = 200


# >>> Helper functions >>>

def _add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        'input_bed_fpath',
        type=str,
        help='Path to input BED file (required)'
    )
    parser.add_argument(
        'input_bam_fpath',
        type=str,
        help='Path to input BAM file (required)'
    )
    parser.add_argument(
        'input_fasta_fpath',
        type=str,
        help='Path to input FASTA file (required)'
    )
    parser.add_argument(
        '--num-read-threshold',
        type=int,
        default=DEFAULT_NUM_READ_THRESHOLD,
        help='Minimum number of read-throughs to consider a region resolved (default: 5)'
    )
    parser.add_argument(
        '--shoulder-len',
        type=int,
        default=DEFAULT_SHOULDER_LEN,
        help='Shoulder length for coordinate checking (default: 200)'
    )
    parser.add_argument(
        '--samtools',
        type=str,
        default=SAMTOOLS_DEFAULT_FPATH,
        help='Path to samtools executable'
    )
    parser.add_argument(
        '--tmpdir',
        type=str,
        default=None,
        help='Temporary directory (default: system temp dir)'
    )
# end def


def _validate_args(args: argparse.Namespace) -> None:

    for fpath in (args.input_bed_fpath,
                  args.input_bam_fpath,
                  args.input_fasta_fpath):
        if not os.path.isfile(fpath):
            sys.stderr.write(f'Error: file `{fpath}` does not exist\n')
            sys.exit(1)
        # end if
    # end for

    if args.num_read_threshold < 1:
        sys.stderr.write(
            'Error: num-read-threshold must be >= 1, got `{}`\n'.format(
                args.num_read_threshold
            )
        )
        sys.exit(1)
    # end if

    if args.shoulder_len < 1:
        sys.stderr.write(
            'Error: shoulder-len must be >= 1, got `{}`\n'.format(
                args.shoulder_len
            )
        )
        sys.exit(1)
    # end if

    _validate_samtools(args.samtools)
# end def


def _validate_samtools(executable_fpath: str) -> None:
    cmd = [
        executable_fpath,
        '--version'
    ]
    pipe = sp.Popen(
        cmd,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        text=True,
        encoding='utf-8',
    )

    _, err = pipe.communicate()
    if pipe.returncode != 0:
        sys.stderr.write(
            'Error: cannot test samtools executable: `{}`\n'.format(
                executable_fpath
            )
        )
        sys.stderr.write('Error: please specify executable with --samtools option\n')
        sys.stderr.write('{}\n'.format(err))
        sys.exit(1)
    # end if
# end def

# <<<


# >>> VerifyArgs class >>>

class VerifyArgs:
    def __init__(self,
                 input_bed_fpath: str,
                 input_bam_fpath: str,
                 input_fasta_fpath: str,
                 num_read_threshold: int = DEFAULT_NUM_READ_THRESHOLD,
                 shoulder_len: int = DEFAULT_SHOULDER_LEN,
                 samtools_fpath: str = SAMTOOLS_DEFAULT_FPATH,
                 tmpdir: Optional[str] = None):
        self.input_bed_fpath: str = input_bed_fpath
        self.input_bam_fpath: str = input_bam_fpath
        self.input_fasta_fpath: str = input_fasta_fpath
        self.num_read_threshold: int = num_read_threshold
        self.shoulder_len: int = shoulder_len
        self.samtools_fpath: str = samtools_fpath
        self.tmpdir: Optional[str] = tmpdir
    # end def

    @classmethod
    def parse_args(cls) -> 'VerifyArgs':
        parser = argparse.ArgumentParser(
            description='Argument parser for reverlor_verify.'
        )
        _add_arguments(parser)
        args = parser.parse_args()
        _validate_args(args)
        # Convert to absolute paths
        args.input_bed_fpath = os.path.abspath(args.input_bed_fpath)
        args.input_bam_fpath = os.path.abspath(args.input_bam_fpath)
        args.input_fasta_fpath = os.path.abspath(args.input_fasta_fpath)
        return cls(
            input_bed_fpath=args.input_bed_fpath,
            input_bam_fpath=args.input_bam_fpath,
            input_fasta_fpath=args.input_fasta_fpath,
            num_read_threshold=args.num_read_threshold,
            shoulder_len=args.shoulder_len,
            samtools_fpath=args.samtools,
            tmpdir=args.tmpdir,
        )
    # end def
# end class

# <<<
