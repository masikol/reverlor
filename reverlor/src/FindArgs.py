#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess as sp
from typing import Optional


MINIMAP2_DEFAULT_FPATH = 'minimap2'
BEDTOOLS_DEFAULT_FPATH = 'bedtools'
MINIMAP_K_DEFAULT = 19
MINIMAP_W_DEFAULT = 19
MINIMAP_M_DEFAULT = 127


# >>> Helper functions >>>

def _add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        'fasta_fpath',
        type=str,
        help='Path to input FASTA file (required)'
    )
    parser.add_argument(
        'output_dir',
        type=str,
        help='Path to output directory (required)'
    )
    parser.add_argument(
        '--min-repeat-len',
        type=int,
        default=200,
        help='Minimum repeat length (default: 200)'
    )
    parser.add_argument(
        '--min-repeat-interval',
        type=int,
        default=100,
        help='Minimum interval between repeats (default: 100). If the interval is shorter, the repeats get merged.'
    )
    parser.add_argument(
        '--minimap2',
        type=str,
        default=MINIMAP2_DEFAULT_FPATH,
        help='Path to minimap2 executable'
    )
    parser.add_argument(
        '--bedtools',
        type=str,
        default=BEDTOOLS_DEFAULT_FPATH,
        help='Path to bedtools executable'
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
# end def


def _validate_args(args: argparse.Namespace) -> None:

    if args.min_repeat_len < 0:
        sys.stderr.write(f'Error: invalid min repeat length: `{args.min_repeat_len}`\n')
        sys.stderr.write('It must be a non-negative integer number\n')
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

    if not os.path.isfile(args.fasta_fpath):
        sys.stderr.write(f'Error: file `{args.fasta_fpath}` does not exist\n')
        sys.exit(1)
    # end if

    _validate_minimap2(args.minimap2)
    _validate_bedtools(args.bedtools)
# end def

def _validate_minimap2(executable_fpath: str) -> None:
    cmd = [
        executable_fpath,
        '--version'
    ]
    pipe = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True, encoding='utf-8')

    _, err = pipe.communicate()
    if pipe.returncode != 0:
        sys.stderr.write(
            'Error: cannot test minimap2 executable: `{}`\n'.format(executable_fpath)
        )
        sys.stderr.write('Error: please specify executable with --minimap2 option\n')
        sys.stderr.write('{}\n'.format(err));
        sys.exit(1)
    # end if
# end def

def _validate_bedtools(executable_fpath: str) -> None:
    cmd = [
        executable_fpath,
        '--version',
    ]
    pipe = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True, encoding='utf-8')

    _, err = pipe.communicate()
    if pipe.returncode != 0:
        sys.stderr.write(
            'Error: cannot test minimap2 executable: `{}`\n'.format(executable_fpath)
        )
        sys.stderr.write('Error: please specify executable with --bedtools option\n')
        sys.stderr.write('{}\n'.format(err));
        sys.exit(1)
    # end if
# end def

# <<<


# >>> FindArgs class >>>

class FindArgs:
    def __init__(self,
                 fasta_fpath: str,
                 output_dir: str,
                 min_repeat_len: int=200,
                 minimap2_fpath: Optional[str]=MINIMAP2_DEFAULT_FPATH,
                 bedtools_fpath: Optional[str]=BEDTOOLS_DEFAULT_FPATH,
                 minimap_k: int=MINIMAP_K_DEFAULT,
                 minimap_w: int=MINIMAP_W_DEFAULT,
                 minimap_m: int=MINIMAP_M_DEFAULT):
        self.fasta_fpath: str = fasta_fpath
        self.output_dir: str = output_dir
        self.min_repeat_len: int = min_repeat_len
        self.minimap2_fpath: Optional[str] = minimap2_fpath
        self.bedtools_fpath: Optional[str] = bedtools_fpath
        self.minimap_k: int = minimap_k
        self.minimap_w: int = minimap_w
        self.minimap_m: int = minimap_m
    # end def

    @classmethod
    def parse_args(cls) -> 'FindArgs':
        parser = argparse.ArgumentParser(
            description='Argument parser for find_repeats.'
        )
        _add_arguments(parser)
        args = parser.parse_args()
        _validate_args(args)
        # Convert to absolute paths
        args.fasta_fpath = os.path.abspath(args.fasta_fpath)
        args.output_dir = os.path.abspath(args.output_dir)
        return cls(
            fasta_fpath=args.fasta_fpath,
            output_dir=args.output_dir,
            min_repeat_len=args.min_repeat_len,
            minimap2_fpath=args.minimap2,
            bedtools_fpath=args.bedtools,
            minimap_k=args.minimap_k,
            minimap_w=args.minimap_w,
            minimap_m=args.minimap_m,
        )
    # end def
# end class

# <<<
