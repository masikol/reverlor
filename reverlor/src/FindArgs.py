#!/usr/bin/env python3

import os
import re
import sys
import argparse
import subprocess as sp
from typing import Optional


MINIMAP2_DEFAULT_FPATH = 'minimap2'
BEDTOOLS_DEFAULT_FPATH = 'bedtools'
MINIMAP_K_DEFAULT = 19
MINIMAP_W_DEFAULT = 19
MINIMAP_M_DEFAULT = 127
MINIMAP_X_CHOICES = ('map-ont', 'lr:hq', 'map-hifi', 'map-pb', 'map-iclr', 'asm5', 'asm10', 'asm20',)
DEFAULT_MIN_REPAT_LEN = 200
DEFAULT_MIN_REPEAT_INTERVAL = 100
FINDER_CHOICES = ('minimap2', 'repeat-scout', 'phraider', 'total-repeats')


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
        default=DEFAULT_MIN_REPAT_LEN,
        help='Minimum repeat length (default: 200)'
    )
    parser.add_argument(
        '--min-repeat-interval',
        type=int,
        default=DEFAULT_MIN_REPEAT_INTERVAL,
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
    parser.add_argument(
        '--minimap-x',
        type=str,
        default=None,
        choices=MINIMAP_X_CHOICES,
        help='minimap2 preset (default: not passed)'
    )
    parser.add_argument(
        '--finder',
        type=str,
        default='minimap2',
        choices=FINDER_CHOICES,
        help='Repeat finder to use (default: minimap2)'
    )
    parser.add_argument(
        '--repeat-scout',
        type=str,
        default=None,
        help='Path to RepeatScout directory containing build_lmer_table and RepeatScout binaries'
    )
    parser.add_argument(
        '--phraider',
        type=str,
        default=None,
        help='Path to phRAIDER executable'
    )
    parser.add_argument(
        '--total-repeats',
        type=str,
        default=None,
        help='Path to TotalRepeats.jar executable'
    )
    parser.add_argument(
        '--java',
        type=str,
        default='java',
        help='Path to Java 25+ executable (default: java)'
    )
    parser.add_argument(
        '--keep-tmp',
        action='store_true',
        default=False,
        help='Keep temporary files (default: False)'
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

    if args.finder == 'repeat-scout':
        _validate_repeat_scout(args.repeat_scout)
    elif args.finder == 'phraider':
        _validate_phraider(args.phraider)
    elif args.finder == 'total-repeats':
        _validate_total_repeats(args.total_repeats)
        _validate_java(args.java)
    else:
        _validate_minimap2(args.minimap2)
    # end if
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
        sys.stderr.write('{}\n'.format(err))
        sys.exit(1)
    # end if
# end def

def _validate_repeat_scout(repeat_scout_dir: Optional[str]) -> None:
    if repeat_scout_dir is None:
        return
    # end if
    if not os.path.isdir(repeat_scout_dir):
        sys.stderr.write(f'Error: RepeatScout directory `{repeat_scout_dir}` does not exist\n')
        sys.exit(1)
    # end if
    for name in ('build_lmer_table', 'RepeatScout'):
        fpath = os.path.join(repeat_scout_dir, name)
        if not os.path.isfile(fpath):
            sys.stderr.write(f'Error: `{name}` not found in `{repeat_scout_dir}`\n')
            sys.exit(1)
        # end if
    # end for
# end def

def _validate_phraider(phraider_fpath: Optional[str]) -> None:
    if phraider_fpath is None:
        return
    # end if
    if not os.path.isfile(phraider_fpath):
        sys.stderr.write(f'Error: phRAIDER executable `{phraider_fpath}` does not exist\n')
        sys.exit(1)
    # end if
    if not os.access(phraider_fpath, os.X_OK):
        sys.stderr.write(f'Error: phRAIDER executable `{phraider_fpath}` is not executable\n')
        sys.exit(1)
    # end if
# end def

def _validate_total_repeats(total_repeats_fpath: Optional[str]) -> None:
    if total_repeats_fpath is None:
        return
    # end if
    if not os.path.isfile(total_repeats_fpath):
        sys.stderr.write(
            'Error: TotalRepeats JAR `{}` does not exist\n'.format(total_repeats_fpath)
        )
        sys.exit(1)
    # end if
# end def

def _validate_java(java_fpath: str) -> None:
    cmd = [java_fpath, '-version']
    pipe = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True, encoding='utf-8')
    out, err = pipe.communicate()
    if pipe.returncode != 0:
        sys.stderr.write(
            'Error: cannot run Java executable `{}`\n'.format(java_fpath)
        )
        sys.stderr.write(err)
        sys.exit(1)
    # end if
    for line in (out or '').splitlines() + (err or '').splitlines():
        m = re.search(r'"(\d+)\.', line)
        if m is not None:
            major = int(m.group(1))
            if major >= 25:
                return
            # end if
            sys.stderr.write(
                'Error: Java version {} is too old, need >= 25\n'.format(major)
            )
            sys.exit(1)
        # end if
    # end for
    sys.stderr.write('Error: could not detect Java version\n')
    sys.exit(1)
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
        sys.stderr.write('{}\n'.format(err))
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
                 minimap_m: int=MINIMAP_M_DEFAULT,
                 minimap_x: Optional[str]=None,
                 finder: str='minimap2',
                 repeat_scout_dir: Optional[str]=None,
                 phraider_fpath: Optional[str]=None,
                 total_repeats_fpath: Optional[str]=None,
                 java_fpath: str='java',
                 keep_tmp: bool=False):
        self.fasta_fpath: str = fasta_fpath
        self.output_dir: str = output_dir
        self.min_repeat_len: int = min_repeat_len
        self.minimap2_fpath: Optional[str] = minimap2_fpath
        self.bedtools_fpath: Optional[str] = bedtools_fpath
        self.minimap_k: int = minimap_k
        self.minimap_w: int = minimap_w
        self.minimap_m: int = minimap_m
        self.minimap_x: Optional[str] = minimap_x
        self.finder: str = finder
        self.repeat_scout_dir: Optional[str] = repeat_scout_dir
        self.phraider_fpath: str = phraider_fpath
        self.total_repeats_fpath: str = total_repeats_fpath
        self.java_fpath: str = java_fpath
        self.keep_tmp: bool = keep_tmp
        if repeat_scout_dir is not None:
            self.build_lmer_table_fpath: str = os.path.join(
                repeat_scout_dir,
                'build_lmer_table'
            )
            self.repeat_scout_fpath: str = os.path.join(
                repeat_scout_dir,
                'RepeatScout'
            )
        else:
            self.build_lmer_table_fpath: str = 'build_lmer_table'
            self.repeat_scout_fpath: str = 'RepeatScout'
        # end if
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
            minimap_x=args.minimap_x,
            finder=args.finder,
            repeat_scout_dir=args.repeat_scout,
            phraider_fpath=args.phraider,
            total_repeats_fpath=args.total_repeats,
            java_fpath=args.java,
            keep_tmp=args.keep_tmp,
        )
    # end def
# end class

# <<<
