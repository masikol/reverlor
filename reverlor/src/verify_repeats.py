#!/usr/bin/env python3

import os
import sys

from Bio import SeqIO

from .VerifyArgs import VerifyArgs
from .CoordIntersecter import CoordIntersecter
from .bed_lib import RepeatRegion, VerifyResult, read_bed_to_regions


def find_unresolved_repeats(args: VerifyArgs) -> list[VerifyResult]:

    os.makedirs(args.output_dir, exist_ok=True)

    all_regions = read_bed_to_regions(args.input_bed_fpath)
    ref_len_dict = _count_ref_lengths(args.input_fasta_fpath)
    coord_intersecter = CoordIntersecter(
        args.input_bam_fpath,
        samtools_fpath=args.samtools_fpath,
        F_flags=[256],
    )

    regions_by_ref = {
        ref_id: [] for ref_id in ref_len_dict
    }
    for region in all_regions:
        regions_by_ref[region.ref_id].append(region)
    # end for

    results: list[VerifyResult] = []

    for ref_id, regions in regions_by_ref.items():
        regions.sort(key=lambda r: r.start)

        for reg_i, region in enumerate(regions):
            shoulder_start = region.start - args.shoulder_len
            shoulder_end   = region.end   + args.shoulder_len

            ref_len = ref_len_dict[region.ref_id]
            skip = False

            if shoulder_start < 1:
                sys.stderr.write(
                    'Region {}: cannot check left shoulder: '
                    'shoulder_start_coord = {} < 1\n'.format(
                        region, shoulder_start
                    )
                )
                skip = True
            # end if
            if shoulder_end > ref_len:
                sys.stderr.write(
                    'Region {}: cannot check right shoulder: '
                    'shoulder_end_coord = {} > ref_len = {}\n'.format(
                        region, shoulder_end, ref_len
                    )
                )
                skip = True
            # end if

            if not skip and _check_coord_within_any_region(
                shoulder_start, regions
            ):
                shoulder_start = _find_shoulder_coord(
                    regions, reg_i, args.shoulder_len, ref_len,
                    right_shoulder=False
                )
                if shoulder_start is None:
                    sys.stderr.write(
                        'Region {}: cannot find left shoulder coord for it: '
                        'they all are within regions or beyond reference\n'
                        .format(region)
                    )
                    skip = True
                # end if
            # end if

            if not skip and _check_coord_within_any_region(
                shoulder_end, regions
            ):
                shoulder_end = _find_shoulder_coord(
                    regions, reg_i, args.shoulder_len, ref_len,
                    right_shoulder=True
                )
                if shoulder_end is None:
                    sys.stderr.write(
                        'Region {}: cannot find right shoulder coord for it: '
                        'they all are within regions or beyond reference\n'
                        .format(region)
                    )
                    skip = True
                # end if
            # end if

            if skip:
                continue
            # end if

            read_ids = coord_intersecter.intersect_coords(
                region.ref_id,
                shoulder_start,
                shoulder_end,
            )

            num_read_throughs = len(read_ids)
            if num_read_throughs < args.num_read_threshold:
                results.append(VerifyResult(
                    region=region,
                    num_read_throughs=num_read_throughs
                ))
                sys.stdout.write(
                    'Region {}: {} read-throughs\n'.format(
                        region,
                        num_read_throughs
                    )
                )
            # end if
        # end for
    # end for

    return results
# end def


def _count_ref_lengths(input_fasta_fpath: str) -> dict[str, int]:
    with open(input_fasta_fpath, 'rt') as ifh:
        seq_records = SeqIO.parse(ifh, 'fasta')
        seq_lengths = {
            sr.id: len(str(sr.seq)) for sr in seq_records
        }
    # end with
    return seq_lengths
# end def


def _check_coord_within_any_region(coord: int,
                                   regions: list[RepeatRegion]) -> bool:
    return any(
        region.start <= coord <= region.end
        for region in regions
    )
# end def


def _find_shoulder_coord(regions: list[RepeatRegion],
                         reg_i: int,
                         shoulder_len: int,
                         ref_len: int,
                         right_shoulder: bool = True) -> int | None:

    if right_shoulder:
        search_regions = regions[reg_i + 1:]
        get_coord = _get_coord_downstream
    else:
        search_regions = list(reversed(regions[:reg_i]))
        get_coord = _get_coord_upstream
    # end if

    for region in search_regions:
        coord = get_coord(region)
        if coord < 0 or coord > ref_len:
            return None
        # end if
        if not _check_coord_within_any_region(coord, search_regions):
            word = 'right' if right_shoulder else 'left'
            sys.stderr.write(
                'INFO: Found {} shoulder coord for region {}: {}\n'.format(
                    word, regions[reg_i], coord
                )
            )
            return coord
        # end if
    # end for

    return None
# end def


def _get_coord_downstream(region: RepeatRegion, shoulder_len: int) -> int:
    return region.end + shoulder_len
# end def


def _get_coord_upstream(region: RepeatRegion, shoulder_len: int) -> int:
    return region.end - shoulder_len
# end def

