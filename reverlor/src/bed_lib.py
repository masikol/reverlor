
import sys
import subprocess as sp

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
            self.start,
            self.end,
            self.end - self.start + 1
        )
    # end def
# end class


def merge_features(args: FindArgs,
                   input_bed_fpath: str,
                   output_bed_fpath: str) -> None:

    sort_cmd = [args.bedtools_fpath, 'sort', '-i', input_bed_fpath]
    merge_cmd = [args.bedtools_fpath, 'merge', '-d', str(args.min_repeat_interval)]

    proc_sort = sp.Popen(
        sort_cmd,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        text=True
    )
    proc_merge = sp.Popen(
        merge_cmd,
        stdin=proc_sort.stdout,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        text=True
    )

    if proc_sort.stdout:
        proc_sort.stdout.close()
    # end if

    out_str, err_str = proc_merge.communicate()
    proc_sort.wait()

    if proc_sort.returncode != 0:
        sys.stderr.write('Error running bedtools sort:\n')
        sys.stderr.write('{}\n'.format(err_str))
        sys.exit(1)
    # end if
    if proc_merge.returncode != 0:
        sys.stderr.write('Error running bedtools merge:\n')
        sys.stderr.write('{}\n'.format(err_str))
        sys.exit(1)
    # end if

    with open(output_bed_fpath, 'w') as out_handle:
        for line in out_str.splitlines():
            values = line.split('\t')
            start_coord = int(values[1]) # 0-based, closed
            end_coord   = int(values[2]) # 0-based, open
            region_len = end_coord - start_coord
            if region_len >= args.min_repeat_len:
                out_handle.write(line)
                out_handle.write('\n')
            # end if
        # end for
    # end with
# end def


def read_bed_to_regions(input_fpath: str) -> list[RepeatRegion]:
    regions = []
    with open(input_fpath, 'rt') as ifh:
        for line in ifh:
            vals = line.strip().split('\t')
            regions.append(RepeatRegion(
                ref_id=vals[0],
                start=int(vals[1]) + 1,   # to 1-based, inclusive
                end=int(vals[2]),         # to 1-based, inclusive
            ))
        # end for
    # end with
    return regions
# end def
