
import sys
import subprocess as sp

from .FindArgs import FindArgs


def merge_features(args: FindArgs,
                   input_bed_fpath: str,
                   output_bed_fpath: str) -> None:

    sort_cmd = [args.bedtools_fpath, 'sort', '-i', input_bed_fpath]
    merge_cmd = [args.bedtools_fpath, 'merge', '-d', '100']

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
