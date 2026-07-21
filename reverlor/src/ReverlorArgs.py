#/usr/bin/env python3

import os
import sys
import argparse
import subprocess as sp
from typing import Optional


MINIMAP2_DEFAULT_FPATH = 'minimap2'
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
        '--minimap2',
        type=str,
        default=MINIMAP2_DEFAULT_FPATH,
        help='Path to minimap2 executable'
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

    _validate_minimap2(args.minimap2)
# end def


def _validate_minimap2(executable_fpath: str) -> None:
    cmd = [executable_fpath, '--version']
    pipe = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True, encoding='utf-8')
    _, err = pipe.communicate()
    if pipe.returncode != 0:
        sys.stderr.write(
            'Error: cannot test minimap2 executable: `{}`\n'.format(executable_fpath)
        )
        sys.stderr.write('Please specify executable with --minimap2 option\n')
        sys.stderr.write('{}\n'.format(err))
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
                 minimap2_fpath: str = MINIMAP2_DEFAULT_FPATH,
                 minimap_k: int = MINIMAP_K_DEFAULT,
                 minimap_w: int = MINIMAP_W_DEFAULT,
                 minimap_m: int = MINIMAP_M_DEFAULT,
                 minimap_x: Optional[str] = None,
                 keep_tmp: bool = False,
                 tmpdir: Optional[str] = None,
                 threads: int = 1):
        self.fasta_fpath: str = fasta_fpath
        self.output_dir: str = output_dir
        self.input_bam_fpath: str = input_bam_fpath
        self.min_repeat_len: int = min_repeat_len
        self.min_repeat_interval: int = min_repeat_interval
        self.num_read_threshold: int = num_read_threshold
        self.shoulder_len: int = shoulder_len
        self.minimap2_fpath: str = minimap2_fpath
        self.minimap_k: int = minimap_k
        self.minimap_w: int = minimap_w
        self.minimap_m: int = minimap_m
        self.minimap_x: Optional[str] = minimap_x
        self.keep_tmp: bool = keep_tmp
        self.tmpdir: Optional[str] = tmpdir
        self.threads: int = threads
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
        return cls(
            fasta_fpath=args.fasta_fpath,
            output_dir=args.output_dir,
            input_bam_fpath=args.input_bam_fpath,
            min_repeat_len=args.min_repeat_len,
            min_repeat_interval=args.min_repeat_interval,
            num_read_threshold=args.num_read_threshold,
            shoulder_len=args.shoulder_len,
            minimap2_fpath=args.minimap2,
            minimap_k=args.minimap_k,
            minimap_w=args.minimap_w,
            minimap_m=args.minimap_m,
            minimap_x=args.minimap_x,
            keep_tmp=args.keep_tmp,
            tmpdir=args.tmpdir,
            threads=args.threads,
        )
    # end def
# end class
