#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess as sp
from typing import Optional

from .ReverlorArgs import ReverlorArgs
from.reverlor_logging import setup_logging


DEFAULT_NUM_READ_THRESHOLD = 5
DEFAULT_SHOULDER_LEN = 200


# >>> Helper functions >>>

def _add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        'input_fasta_fpath',
        type=str,
        help='Path to input FASTA file (required)'
    )
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
        'output_dir',
        type=str,
        help='Path to output directory (required)'
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
        '--tmpdir',
        type=str,
        default=None,
        help='Temporary directory (default: system temp dir)'
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action='count',
        default=0,
        help='Increase verbosity. Can be used multiple times: -v, -vv, up to -vvv, which is debug mode.'
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
# end def

# <<<


# >>> VerifyArgs class >>>

class VerifyArgs:
    def __init__(self,
                 input_fasta_fpath: str,
                 input_bed_fpath: str,
                 input_bam_fpath: str,
                 output_dir: str,
                 num_read_threshold: int = DEFAULT_NUM_READ_THRESHOLD,
                 shoulder_len: int = DEFAULT_SHOULDER_LEN,
                 tmpdir: Optional[str] = None):
        self.input_fasta_fpath: str = input_fasta_fpath
        self.input_bed_fpath: str = input_bed_fpath
        self.input_bam_fpath: str = input_bam_fpath
        self.output_dir: str = output_dir
        self.num_read_threshold: int = num_read_threshold
        self.shoulder_len: int = shoulder_len
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
        args.input_fasta_fpath = os.path.abspath(args.input_fasta_fpath)
        args.input_bed_fpath = os.path.abspath(args.input_bed_fpath)
        args.input_bam_fpath = os.path.abspath(args.input_bam_fpath)
        args.output_dir = os.path.abspath(args.output_dir)
        setup_logging(args.verbose)
        return cls(
            input_fasta_fpath=args.input_fasta_fpath,
            input_bed_fpath=args.input_bed_fpath,
            input_bam_fpath=args.input_bam_fpath,
            output_dir=args.output_dir,
            num_read_threshold=args.num_read_threshold,
            shoulder_len=args.shoulder_len,
            tmpdir=args.tmpdir,
        )
    # end def

    @classmethod
    def from_reverlor_args(cls,
                           rev: 'ReverlorArgs',
                           input_bed_fpath: str) -> 'VerifyArgs':
        if not os.path.isfile(input_bed_fpath):
            sys.stderr.write(
                'Error: input BED file `{}` does not exist\n'.format(
                    input_bed_fpath
                )
            )
            sys.exit(1)
        # end if
        return cls(
            input_fasta_fpath=rev.fasta_fpath,
            input_bed_fpath=input_bed_fpath,
            input_bam_fpath=rev.input_bam_fpath,
            output_dir=rev.output_dir,
            num_read_threshold=rev.num_read_threshold,
            shoulder_len=rev.shoulder_len,
            tmpdir=rev.tmpdir,
        )
    # end def
# end class

# <<<
