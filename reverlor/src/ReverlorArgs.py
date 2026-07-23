#/usr/bin/env python3

import os
import sys
import argparse
import subprocess as sp
from typing import Optional

import src.defaults as defaults
from .reverlor_logging import setup_logging
from ._version import __version__, __last_update_date__


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
        'fasta_fpath',
        type=str,
        help='Path to input FASTA file'
    )
    parser.add_argument(
        'input_bam_fpath',
        type=str,
        help='Path to input BAM file for verification'
    )
    parser.add_argument(
        'output_dir',
        type=str,
        help='Path to output directory'
    )
    parser.add_argument(
        '-l',
        '--min-repeat-len',
        type=int,
        default=defaults.MIN_REPAT_LEN,
        help=f'Minimum repeat length (default: {defaults.MIN_REPAT_LEN})'
    )
    parser.add_argument(
        '-i',
        '--min-repeat-interval',
        type=int,
        default=defaults.MIN_REPEAT_INTERVAL,
        help=(
            f'Minimum interval between repeats (default: {defaults.MIN_REPEAT_INTERVAL}).'
            'If the interval is shorter, the repeats get merged.'
        )
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
        '-k',
        '--minimap-k',
        type=int,
        default=defaults.MINIMAP_K,
        help=f'minimap2 k-mer length (default: {defaults.MINIMAP_K})'
    )
    parser.add_argument(
        '-w',
        '--minimap-w',
        type=int,
        default=defaults.MINIMAP_W,
        help=f'minimap2 minimizer window size (default: {defaults.MINIMAP_W})'
    )
    parser.add_argument(
        '-m',
        '--minimap-m',
        type=int,
        default=defaults.MINIMAP_M,
        help=f'minimap2 matching score (default: {defaults.MINIMAP_M})'
    )
    parser.add_argument(
        '--minimap-x',
        type=str,
        default=defaults.MINIMAP_X,
        choices=defaults.MINIMAP_X_CHOICES,
        help='minimap2 preset (default: not set)'
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
        '--keep-tmp',
        action='store_true',
        default=defaults.KEEP_TMP_FILES,
        help='Keep temporary files (default: {defaults.KEEP_TMP_FILES})'
    )
    parser.add_argument(
        '--tmpdir',
        type=str,
        default=defaults.TMP_DIR_PATH,
        help='Temporary directory (default: system temp dir)'
    )
    parser.add_argument(
        '-t',
        '--threads',
        type=int,
        default=defaults.NUM_THREADS,
        help=f'number of CPU threads to use (default: {defaults.NUM_THREADS})'
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

    if not os.path.isfile(args.fasta_fpath):
        sys.stderr.write(f'Error: file `{args.fasta_fpath}` does not exist\n')
        sys.exit(1)
    # end if

    if not os.path.isfile(args.input_bam_fpath):
        sys.stderr.write(f'Error: file `{args.input_bam_fpath}` does not exist\n')
        sys.exit(1)
    # end if

    if args.min_repeat_len < 0:
        sys.stderr.write(f'Error: invalid min repeat length: `{args.min_repeat_len}`\n')
        sys.stderr.write('It must be a non-negative integer number\n')
        sys.exit(1)
    # end if

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

    if args.minimap_k < 0:
        sys.stderr.write(f'Error: invalid minimap-k: `{args.minimap_k}`\n')
        sys.stderr.write('It must be a non-negative integer\n')
        sys.exit(1)
    # end if

    if args.minimap_w < 0:
        sys.stderr.write(f'Error: invalid minimap-w: `{args.minimap_w}`\n')
        sys.stderr.write('It must be a non-negative integer\n')
        sys.exit(1)
    # end if

    if args.minimap_m <= 0:
        sys.stderr.write(f'Error: invalid minimap-m: `{args.minimap_m}`\n')
        sys.stderr.write('It must be a positive integer\n')
        sys.exit(1)
    # end if
# end def


class ReverlorArgs:
    def __init__(self,
                 fasta_fpath: str,
                 output_dir: str,
                 input_bam_fpath: str,
                 min_repeat_len: int = defaults.MIN_REPAT_LEN,
                 min_repeat_interval: int = defaults.MIN_REPEAT_INTERVAL,
                 span_threshold: int = defaults.NUM_READ_THRESHOLD,
                 shoulder_len: int = defaults.SHOULDER_LEN,
                 minimap_k: int = defaults.MINIMAP_K,
                 minimap_w: int = defaults.MINIMAP_W,
                 minimap_m: int = defaults.MINIMAP_M,
                 minimap_x: Optional[str] = None,
                 keep_tmp: bool = False,
                 samtools_f: Optional[list[int]] = None,
                 samtools_F: Optional[list[int]] = None,
                 tmpdir: Optional[str] = None,
                 threads: int = 1):
        self.fasta_fpath: str = fasta_fpath
        self.output_dir: str = output_dir
        self.input_bam_fpath: str = input_bam_fpath
        self.min_repeat_len: int = min_repeat_len
        self.min_repeat_interval: int = min_repeat_interval
        self.span_threshold: int = span_threshold
        self.shoulder_len: int = shoulder_len
        self.minimap_k: int = minimap_k
        self.minimap_w: int = minimap_w
        self.minimap_m: int = minimap_m
        self.minimap_x: Optional[str] = minimap_x
        self.keep_tmp: bool = keep_tmp
        self.samtools_f: list[int] = samtools_f if samtools_f is not None else defaults.SAMTOOLS_f
        self.samtools_F: list[int] = samtools_F if samtools_F is not None else defaults.SAMTOOLS_F
        self.tmpdir: Optional[str] = tmpdir
        self.threads: int = threads
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
    def parse_args(cls) -> 'ReverlorArgs':
        parser = argparse.ArgumentParser(
            description='Reverlor: repeat verification using long reads.'
        )
        _add_arguments(parser)
        args = parser.parse_args()
        _validate_args(args)
        args.fasta_fpath = os.path.abspath(args.fasta_fpath)
        args.output_dir = os.path.abspath(args.output_dir)
        args.input_bam_fpath = os.path.abspath(args.input_bam_fpath)
        setup_logging(args.verbose)
        return cls(
            fasta_fpath=args.fasta_fpath,
            output_dir=args.output_dir,
            input_bam_fpath=args.input_bam_fpath,
            min_repeat_len=args.min_repeat_len,
            min_repeat_interval=args.min_repeat_interval,
            span_threshold=args.span,
            shoulder_len=args.shoulder_len,
            minimap_k=args.minimap_k,
            minimap_w=args.minimap_w,
            minimap_m=args.minimap_m,
            minimap_x=args.minimap_x,
            keep_tmp=args.keep_tmp,
            samtools_f=args.samtools_f,
            samtools_F=args.samtools_F,
            tmpdir=args.tmpdir,
            threads=args.threads,
        )
    # end def
# end class
