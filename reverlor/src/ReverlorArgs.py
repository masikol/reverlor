#/usr/bin/env python3

import os
import sys
import argparse
import subprocess as sp
from typing import Optional

from ._version import __version__, __last_update_date__

from.reverlor_logging import setup_logging


MINIMAP_K_DEFAULT = 19
MINIMAP_W_DEFAULT = 19
MINIMAP_M_DEFAULT = 65
MINIMAP_X_CHOICES = ('map-ont', 'lr:hq', 'map-hifi', 'map-pb', 'map-iclr', 'asm5', 'asm10', 'asm20',)
DEFAULT_MIN_REPAT_LEN = 200
DEFAULT_MIN_REPEAT_INTERVAL = 100
DEFAULT_NUM_READ_THRESHOLD = 5
DEFAULT_SHOULDER_LEN = 200


def _add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '-V', '--version',
        action='version',
        version='%(prog)s ' + __version__ + (', ' + __last_update_date__ + ' edition' if __last_update_date__ else ''),
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
        '--min-repeat-len',
        type=int,
        default=DEFAULT_MIN_REPAT_LEN,
        help='Minimum repeat length (default: 200)'
    )
    parser.add_argument(
        '--min-repeat-interval',
        type=int,
        default=DEFAULT_MIN_REPEAT_INTERVAL,
        help='Minimum interval between repeats (default: 100)'
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
        '--minimap-k',
        type=int,
        default=MINIMAP_K_DEFAULT,
        help='minimap2 k-mer length (default: 19)'
    )
    parser.add_argument(
        '--minimap-w',
        type=int,
        default=MINIMAP_W_DEFAULT,
        help='minimap2 minimizer window size (default: 19)'
    )
    parser.add_argument(
        '--minimap-m',
        type=int,
        default=MINIMAP_M_DEFAULT,
        help='minimap2 matching score (default: 127)'
    )
    parser.add_argument(
        '--minimap-x',
        type=str,
        default=None,
        choices=MINIMAP_X_CHOICES,
        help='minimap2 preset (default: not passed)'
    )
    parser.add_argument(
        '--samtools-f',
        type=int,
        action='append',
        default=None,
        help='samtools -f flag (repeatable)'
    )
    parser.add_argument(
        '--samtools-F',
        type=int,
        action='append',
        default=None,
        help='samtools -F flag (repeatable, default: 256)'
    )
    parser.add_argument(
        '--keep-tmp',
        action='store_true',
        default=False,
        help='Keep temporary files (default: False)'
    )
    parser.add_argument(
        '--tmpdir',
        type=str,
        default=None,
        help='Temporary directory (default: system temp dir)'
    )
    parser.add_argument(
        '-t',
        '--threads',
        type=int,
        default=1,
        help='number of CPU threads to use'
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
                 min_repeat_len: int = DEFAULT_MIN_REPAT_LEN,
                 min_repeat_interval: int = DEFAULT_MIN_REPEAT_INTERVAL,
                 num_read_threshold: int = DEFAULT_NUM_READ_THRESHOLD,
                 shoulder_len: int = DEFAULT_SHOULDER_LEN,
                 minimap_k: int = MINIMAP_K_DEFAULT,
                 minimap_w: int = MINIMAP_W_DEFAULT,
                 minimap_m: int = MINIMAP_M_DEFAULT,
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
        self.num_read_threshold: int = num_read_threshold
        self.shoulder_len: int = shoulder_len
        self.minimap_k: int = minimap_k
        self.minimap_w: int = minimap_w
        self.minimap_m: int = minimap_m
        self.minimap_x: Optional[str] = minimap_x
        self.keep_tmp: bool = keep_tmp
        self.samtools_f: list[int] = samtools_f if samtools_f is not None else []
        self.samtools_F: list[int] = samtools_F if samtools_F is not None else [256,]
        self.tmpdir: Optional[str] = tmpdir
        self.threads: int = threads
    # end def

    def __str__(self) -> str:
        lines = ['ReverlorArgs:']
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
            num_read_threshold=args.num_read_threshold,
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
