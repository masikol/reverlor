#!/usr/bin/env python3

import os
import sys
import subprocess as sp
from .FindArgs import FindArgs



def find_repeats(args: FindArgs) -> str:

    input_fasta = os.path.abspath(args.fasta_fpath)
    output_dir = os.path.abspath(args.output_dir)

    os.makedirs(output_dir, exist_ok=True)

    # Process no_rotate mode
    no_rotate_dir = os.path.join(output_dir, 'no_rotate')
    no_rotate_bed = _find_repeats_single(
        input_fasta=input_fasta,
        window_len=args.window_size,
        min_repeat_len=args.min_repeat_len,
        rotate=False,
        output_dir=no_rotate_dir,
        args=args
    )

    # Process rotate mode
    rotate_dir = os.path.join(output_dir, 'rotate')
    rotate_bed = _find_repeats_single(
        input_fasta=input_fasta,
        window_len=args.window_size,
        min_repeat_len=args.min_repeat_len,
        rotate=True,
        output_dir=rotate_dir,
        args=args
    )

    # Merge BED files
    merged_bed = _merge_bed_files(
        no_rotate_bed,
        rotate_bed,
        output_dir,
        args.bedtools_fpath
    )

    print('', file=sys.stderr)
    print('INFO: Completed!', file=sys.stderr)
    print(f'INFO: output directory: `{output_dir}`', file=sys.stderr)
    print(f'INFO: main output file: `{merged_bed}`', file=sys.stderr)

    return merged_bed
# end def



# >>> Helper functions >>>

def _find_repeats_single(input_fasta: str,
                         window_len: int,
                         min_repeat_len: int,
                         rotate: bool,
                         output_dir: str,
                         args: FindArgs) -> str:
    os.makedirs(output_dir, exist_ok=True)

    sliding_chunks_fasta = os.path.join(output_dir, 'sliding_chunks.fasta.gz')
    srt_bam = os.path.join(output_dir, 'repeats.srt.bam')
    repeats_bed = os.path.join(output_dir, 'repeats.bed')

    _make_sliding_chunks(
        input_fasta,
        window_len,
        rotate,
        sliding_chunks_fasta,
        args.seqkit_fpath
    )
    _map_sliding_chunks(
        input_fasta,
        sliding_chunks_fasta,
        srt_bam,
        args.minimap2_fpath,
        args.samtools_fpath
    )
    _make_repeats_bed(
        input_fasta,
        srt_bam,
        min_repeat_len,
        repeats_bed,
        args.con_hi_fpath
    )

    return repeats_bed
# end def


def _make_sliding_chunks(input_fasta: str,
                         window_len: int,
                         rotate: bool,
                         output_fasta: str,
                         seqkit_fpath: str) -> None:

    if rotate:
        half_window_len = window_len // 2
        cmd1 = [
            seqkit_fpath, 'restart',
            '-i', str(half_window_len),
            input_fasta
        ]
        cmd2 = [
            seqkit_fpath, 'sliding',
            '-j', '1',
            '-s', str(window_len),
            '-W', str(window_len)
        ]
        p1 = sp.Popen(cmd1, stdout=sp.PIPE)
        p2 = sp.Popen(cmd2, stdin=p1.stdout, stdout=sp.PIPE)
        p1.stdout.close()
    else:
        cmd = [
            seqkit_fpath,
            'sliding',
            '-j', '1',
            '-s', str(window_len),
            '-W', str(window_len),
            input_fasta
        ]
        p2 = sp.Popen(cmd, stdout=sp.PIPE)
    # end if

    with open(output_fasta, 'wb') as f:
        sp.run(['gzip'], stdin=p2.stdout, stdout=f)
    # end with

    p2.wait()
    if rotate:
        p1.wait()
    # end if
# end def


def _run_command(cmd: list, **kwargs) -> None:
    sp.run(cmd, check=True, **kwargs)
# end def


def _map_sliding_chunks(input_fasta: str,
                        sliding_chunks_fasta: str,
                        srt_bam: str,
                        minimap2_fpath: str,
                        samtools_fpath: str) -> None:

    cmd1 = [
        minimap2_fpath,
        '-a',
        '-x', 'asm20',
        '--eqx',
        '-t', '1',
        input_fasta,
        sliding_chunks_fasta
    ]
    cmd2 = [samtools_fpath, 'view', '-O', 'BAM']
    cmd3 = [samtools_fpath, 'sort', '-O', 'BAM', '-o', srt_bam]

    p1 = sp.Popen(cmd1, stdout=sp.PIPE)
    p2 = sp.Popen(cmd2, stdin=p1.stdout, stdout=sp.PIPE)
    p1.stdout.close()
    p3 = sp.Popen(cmd3, stdin=p2.stdout)
    p2.stdout.close()

    p3.wait()
    p2.wait()
    p1.wait()

    cmd4 = [samtools_fpath, 'index', srt_bam]
    _run_command(cmd4)
# end def


def _make_repeats_bed(input_fasta: str,
                      srt_bam: str,
                      min_repeat_len: int,
                      output_bed: str,
                      con_hi_fpath: str) -> None:

    cmd = [
        'python3',
        con_hi_fpath,
        '-f', input_fasta,
        '-b', srt_bam,
        '-o', output_bed,
        '-O', 'bed',
        '-c', 'off',
        '--no-zero-output',
        '-X', 'off',
        '-C', '1',
        '--min-feature-len', str(min_repeat_len)
    ]
    _run_command(cmd)
# end def


def _merge_bed_files(no_rotate_bed: str,
                     rotate_bed: str,
                     output_dir: str,
                     bedtools_fpath: str) -> str:
    concat_bed = os.path.join(output_dir, 'repeats_concat.bed')
    concat_srt_bed = os.path.join(output_dir, 'repeats_concat.srt.bed')
    merged_bed = os.path.join(output_dir, 'repeats_final.bed')

    with open(concat_bed, 'w') as out_handle:
        for bed_file in [no_rotate_bed, rotate_bed]:
            with open(bed_file, 'r') as in_handle:
                for line in in_handle:
                    out_handle.write(line)
                # end for
            # end with
        # end for
    # end with

    # Sort and merge
    sort_cmd = [bedtools_fpath, 'sort', '-i', concat_bed]
    with open(concat_srt_bed, 'w') as out_handle:
        _run_command(sort_cmd, stdout=out_handle)
    # enw with

    merge_cmd = [bedtools_fpath, 'merge', '-i', concat_srt_bed]

    with open(merged_bed, 'w') as out_handle:
        _run_command(merge_cmd, stdout=out_handle)
    # end with

    # Clean up temporary files
    for tmp_file in [concat_bed, concat_srt_bed]:
        if os.path.isfile(tmp_file):
            os.remove(tmp_file)
        # end if
    # end for

    return merged_bed
# end def

# <<<
