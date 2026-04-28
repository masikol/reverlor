#!/usr/bin/env python3

import sys
import os
import argparse
from typing import Optional


SEQKIT_DEFAULT_FPATH = 'seqkit'
MINIMAP2_DEFAULT_FPATH = 'minimap2'
SAMTOOLS_DEFAULT_FPATH = 'samtools'
CON_HI_DEFAULT_FPATH = 'con-hi.py'
BEDTOOLS_DEFAULT_FPATH = 'bedtools'


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
        '--window-size',
        type=int,
        default=1000,
        help='Window size (default: 1000)'
    )
    parser.add_argument(
        '--min-repeat-len',
        type=int,
        default=200,
        help='Minimum repeat length (default: 200)'
    )
    parser.add_argument(
        '--seqkit',
        type=str,
        default=SEQKIT_DEFAULT_FPATH,
        help='Path to seqkit executable'
    )
    parser.add_argument(
        '--minimap2',
        type=str,
        default=MINIMAP2_DEFAULT_FPATH,
        help='Path to minimap2 executable'
    )
    parser.add_argument(
        '--samtools',
        type=str,
        default=SAMTOOLS_DEFAULT_FPATH,
        help='Path to samtools executable'
    )
    parser.add_argument(
        '--con-hi',
        type=str,
        default=CON_HI_DEFAULT_FPATH,
        help='Path to con-hi.py executable'
    )
    parser.add_argument(
        '--bedtools',
        type=str,
        default=BEDTOOLS_DEFAULT_FPATH,
        help='Path to bedtools executable'
    )
# end def


def _validate_args(args: argparse.Namespace) -> None:
    # Check window_size >= 0
    if args.window_size < 0:
        print(f'Error: invalid window size: `{args.window_size}`', file=sys.stderr)
        print('It must be a non-negative integer number', file=sys.stderr)
        sys.exit(1)
    # end if

    # Check min_repeat_len >= 0
    if args.min_repeat_len < 0:
        print(f'Error: invalid min repeat length: `{args.min_repeat_len}`', file=sys.stderr)
        print('It must be a non-negative integer number', file=sys.stderr)
        sys.exit(1)
    # end if

    # Check input FASTA file exists
    if not os.path.isfile(args.fasta_fpath):
        print(f'Error: file `{args.fasta_fpath}` does not exist', file=sys.stderr)
        sys.exit(1)
    # end if
# end def

# <<<


# >>> FindArgs class >>>

class FindArgs:
    def __init__(self,
                 fasta_fpath: str,
                 output_dir: str,
                 window_size: int=1000,
                 min_repeat_len: int=200,
                 seqkit_fpath: Optional[str]=SEQKIT_DEFAULT_FPATH,
                 minimap2_fpath: Optional[str]=MINIMAP2_DEFAULT_FPATH,
                 samtools_fpath: Optional[str]=SAMTOOLS_DEFAULT_FPATH,
                 con_hi_fpath: Optional[str]=CON_HI_DEFAULT_FPATH,
                 bedtools_fpath: Optional[str]=BEDTOOLS_DEFAULT_FPATH):
        self.fasta_fpath: str = fasta_fpath
        self.output_dir: str = output_dir
        self.window_size: int = window_size
        self.min_repeat_len: int = min_repeat_len
        self.seqkit_fpath: Optional[str] = seqkit_fpath
        self.minimap2_fpath: Optional[str] = minimap2_fpath
        self.samtools_fpath: Optional[str] = samtools_fpath
        self.con_hi_fpath: Optional[str] = con_hi_fpath
        self.bedtools_fpath: Optional[str] = bedtools_fpath
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
            window_size=args.window_size,
            min_repeat_len=args.min_repeat_len,
            seqkit_fpath=args.seqkit,
            minimap2_fpath=args.minimap2,
            samtools_fpath=args.samtools,
            con_hi_fpath=args.con_hi,
            bedtools_fpath=args.bedtools
        )
    # end def
# end class

# <<<
