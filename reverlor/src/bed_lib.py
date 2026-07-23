
from typing import NamedTuple

import pybedtools

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
    (
        pybedtools.BedTool(input_bed_fpath)
            .sort()
            .merge(d=args.min_repeat_interval)
            .filter(lambda ivl: (ivl.end - ivl.start) >= args.min_repeat_len)
            .saveas(output_bed_fpath)
    )
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
