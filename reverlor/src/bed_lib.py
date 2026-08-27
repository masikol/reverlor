
from typing import NamedTuple

from .FindArgs import FindArgs


class RepeatRegion:
    def __init__(self, ref_id, start, end):
        self.ref_id = ref_id
        self.start  = start
        self.end    = end
    # end def
    def __str__(self):
        return '{}:{}-{} (len {:,})'.format(
            self.ref_id,
            self.start+1, # to 1-based, closed
            self.end,     # to 1-based, closed
            self.end - self.start
        )
    # end def
# end class


class VerifyResult(NamedTuple):
    region: RepeatRegion
    num_read_throughs: int
# end class


def merge_features(args: FindArgs,
                   input_bed_fpath: str,
                   output_bed_fpath: str) -> None:
    regions = read_bed_to_regions(input_bed_fpath)
    regions = sort_regions(regions)
    regions = merge_regions(regions, args.min_repeat_interval)
    regions = filter_regions(regions, args.min_repeat_len)
    with open(output_bed_fpath, 'w') as fh:
        for region in regions:
            fh.write('\t'.join((
                region.ref_id,
                str(region.start),
                str(region.end),
            )) + '\n')
        # end for
    # end with
# end def


def sort_regions(regions: list[RepeatRegion]) -> list[RepeatRegion]:
    return sorted(regions, key=lambda r: (r.ref_id, r.start, r.end))
# end def


def merge_regions(regions: list[RepeatRegion],
                  min_repeat_interval: int) -> list[RepeatRegion]:
    merged: list[RepeatRegion] = []
    for region in regions:
        merge_to_prev = len(merged) != 0 \
                        and merged[-1].ref_id == region.ref_id \
                        and region.start - merged[-1].end <= min_repeat_interval
        if merge_to_prev:
            prev = merged[-1]
            prev.end = max(prev.end, region.end)
        else:
            merged.append(RepeatRegion(
                ref_id=region.ref_id,
                start=region.start,
                end=region.end,
            ))
        # end if
    # end for
    return merged
# end def


def filter_regions(regions: list[RepeatRegion],
                   min_repeat_len: int) -> list[RepeatRegion]:
    return [r for r in regions if (r.end - r.start) >= min_repeat_len]
# end def


def read_bed_to_regions(input_fpath: str) -> list[RepeatRegion]:
    regions = []
    with open(input_fpath, 'rt') as ifh:
        for line in ifh:
            vals = line.strip().split('\t')
            regions.append(RepeatRegion(
                ref_id=vals[0].strip(),
                start=int(vals[1].strip()),   # keep 0-based, close
                end=int(vals[2].strip()),     # keep 0-based, open
            ))
        # end for
    # end with
    return regions
# end def


def verify_results_to_bed(verify_results: list[VerifyResult],
                          out_fpath: str) -> None:
    with open(out_fpath, 'w') as fh:
        for vr in verify_results:
            fh.write('\t'.join((
                vr.region.ref_id,
                str(vr.region.start),
                str(vr.region.end),
                'repeat',
                str(vr.num_read_throughs),
            )) + '\n')
        # end for
    # end with
# end def
