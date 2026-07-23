#!/usr/bin/env python3

import os
import sys

import mappy as mp

from .FindArgs import FindArgs
from .bed_lib import merge_features
from .util import rm_files_if_exist


def find_repeats(args: FindArgs) -> str:

    os.makedirs(args.output_dir, exist_ok=True)

    raw_bed_fpath = os.path.join(args.output_dir, 'repeats_raw.bed')
    merged_bed_fpath = os.path.join(args.output_dir, 'repeats.bed')

    _create_raw_repeat_file(args, raw_bed_fpath)
    merge_features(args, raw_bed_fpath, merged_bed_fpath)

    if not args.keep_tmp:
        rm_files_if_exist(raw_bed_fpath)
    # end if

    return merged_bed_fpath
# end def


def _create_raw_repeat_file(args: FindArgs,
                            output_bed_fpath: str) -> None:
    # See files main.c and minimap.h of minimap2
    MM_F_ALL_CHAINS = 0x800000 # for -P
    MM_F_NO_DIAG    =    0x001 # for -D 
    MM_OUT_CG       =    0x020 # for -c, won’t take effect, mappy enforces it anyway, but just in case
    MM_F_CIGAR      =    0x004 # for -c, won’t take effect, mappy enforces it anyway, but just in case

    extra_flags = MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_OUT_CG | MM_F_CIGAR

    aligner = mp.Aligner(
        args.fasta_fpath,
        preset=args.minimap_x,
        k=args.minimap_k,
        w=args.minimap_w,
        min_chain_score=args.minimap_m,
        n_threads=args.threads,
        extra_flags=extra_flags
    )

    with open(output_bed_fpath, 'wt') as bed_handle:
        for name, seq, qual in mp.fastx_read(args.fasta_fpath):
            # Passing name to aligner.map is neccessary for MM_F_NO_DIAG to actually take affect
            for hit in aligner.map(seq, name=name):
                bed_handle.write(
                    _make_bed_string(hit, name)
                )
            # end for
        # end for
    # end with
# end def

def _make_bed_string(hit: mp.Alignment, query_name: str) -> str:
    return '{}\n'.format('\t'.join([
        query_name,
        str(hit.q_st),
        str(hit.q_en),
        hit.ctg,
        str(hit.r_st),
        str(hit.r_en),
    ]))
# end def
