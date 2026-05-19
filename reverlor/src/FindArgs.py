#!/usr/bin/env python3

import sys
import os
import argparse
from typing import Optional


MINIMAP2_DEFAULT_FPATH = 'minimap2'
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
# end def


def _validate_args(args: argparse.Namespace) -> None:

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
                 min_repeat_len: int=200,
                 minimap2_fpath: Optional[str]=MINIMAP2_DEFAULT_FPATH,
                 bedtools_fpath: Optional[str]=BEDTOOLS_DEFAULT_FPATH):
        self.fasta_fpath: str = fasta_fpath
        self.output_dir: str = output_dir
        self.min_repeat_len: int = min_repeat_len
        self.minimap2_fpath: Optional[str] = minimap2_fpath
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
            min_repeat_len=args.min_repeat_len,
            minimap2_fpath=args.minimap2,
            bedtools_fpath=args.bedtools
        )
    # end def
# end class

# <<<
