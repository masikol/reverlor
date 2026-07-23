#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess as sp
from typing import Optional

import src.defaults as defaults
from .ReverlorArgs import ReverlorArgs
from .reverlor_logging import setup_logging
from ._version import __version__, __last_update_date__




# >>> Helper functions >>>

def _add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '-V', '--version',
        action='version',
        version=(
            '%(prog)s ' + __version__ +
            ', ' + __last_update_date__ + ' edition'
        ),
    )
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
        '-s',
        '--span',
        type=int,
        default=defaults.NUM_READ_THRESHOLD,
        help=f'Minimum number of spanning reads to consider a repeat resolved (default: {defaults.NUM_READ_THRESHOLD})'
    )
    parser.add_argument(
        '-u',
        '--shoulder-len',
        type=int,
        default=defaults.SHOULDER_LEN,
        help=f'Shoulder length for coordinate checking (default: {defaults.SHOULDER_LEN})'
    )
    parser.add_argument(
        '--samtools-f',
        type=int,
        action='append',
        default=None,
        help='samtools -f flag (repeatable, default: none set)'
    )
    parser.add_argument(
        '--samtools-F',
        type=int,
        action='append',
        default=None,
        help=f'samtools -F flag (repeatable, default: {defaults.SAMTOOLS_F})'
    )
    parser.add_argument(
        '--tmpdir',
        type=str,
        default=defaults.TMP_DIR_PATH,
        help='Temporary directory (default: system temp dir)'
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action='count',
        default=defaults.VERBOSE,
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

    if args.span < 1:
        sys.stderr.write(
            'Error: --span must be >= 1, got `{}`\n'.format(
                args.span
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
                 span_threshold: int = defaults.NUM_READ_THRESHOLD,
                 shoulder_len: int = defaults.SHOULDER_LEN,
                 samtools_f: Optional[list[int]] = None,
                 samtools_F: Optional[list[int]] = None,
                 tmpdir: Optional[str] = None):
        self.input_fasta_fpath: str = input_fasta_fpath
        self.input_bed_fpath: str = input_bed_fpath
        self.input_bam_fpath: str = input_bam_fpath
        self.output_dir: str = output_dir
        self.span_threshold: int = span_threshold
        self.shoulder_len: int = shoulder_len
        self.samtools_f: list[int] = samtools_f if samtools_f is not None else defaults.SAMTOOLS_f
        self.samtools_F: list[int] = samtools_F if samtools_F is not None else defaults.SAMTOOLS_F
        self.tmpdir: Optional[str] = tmpdir
    # end def

    def __str__(self) -> str:
        lines = ['Run parameters:']
        for attr, val in self.__dict__.items():
            lines.append(f'  {attr:25s}= {val}')
        # end for
        lines.append('='*20)
        return '\n'.join(lines)
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
            span_threshold=args.span,
            shoulder_len=args.shoulder_len,
            samtools_f=args.samtools_f,
            samtools_F=args.samtools_F,
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
            span_threshold=rev.span_threshold,
            shoulder_len=rev.shoulder_len,
            samtools_f=rev.samtools_f,
            samtools_F=rev.samtools_F,
            tmpdir=rev.tmpdir,
        )
    # end def
# end class

# <<<
