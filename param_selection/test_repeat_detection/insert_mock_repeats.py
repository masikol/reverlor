#!/usr/bin/env python3

import re
import os
import sys
import random

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from mock_repeats_settings import RATE_FROM, \
                                  RATE_TO, \
                                  RATE_STEP, \
                                  RANDOM_SEED


random.seed(RANDOM_SEED)

REPEAT_COORDS_RE = re.compile(r'start=([1-9][0-9]*) end=([1-9][0-9]*)')
RECORD_ID_RE = re.compile(r'MR_l(\d+)_r(\d+)_c(\d+)_t(\S+)')

IN_GENOME_FPATH = os.path.abspath(sys.argv[1])
IN_PLASMID_FPATH = os.path.abspath(sys.argv[2])
IN_MOCK_REPEATS_DIRPATH = os.path.abspath(sys.argv[3])
N_DESIRED_REPEAT_COPIES = int(sys.argv[4])

assert N_DESIRED_REPEAT_COPIES > 1
N_REPEAT_COPIES_TO_INSERT = N_DESIRED_REPEAT_COPIES - 1

IN_MOCK_REPEATS_FPATH = os.path.join(
    IN_MOCK_REPEATS_DIRPATH,
    'mock_repeats.fasta'
)

SHOULDER_LEN = 500



def make_mock_repat_fasta_fpath(rate=None):
    if rate is None:
        return IN_MOCK_REPEATS_FPATH
    # end if

    rate_str = '{:.2f}'.format(rate).replace('.', 'dot')
    out_basename = os.path.basename(IN_MOCK_REPEATS_FPATH).replace(
        '.fasta',
        '_{}.fasta'.format(rate_str)
    )
    return os.path.join(
        IN_MOCK_REPEATS_DIRPATH,
        out_basename
    )
# end def

def make_chr_with_insert_base_name(mock_repeat_seq_record, rate=None):
    basic_record_id = str(mock_repeat_seq_record.id)

    rate_str = '{:.2f}'.format(0)
    if rate is not None:
        rate_str = '{:.2f}'.format(rate)
    # end if

    match = re.match(RECORD_ID_RE, basic_record_id)
    l, r, c, t = match.group(1), match.group(2), match.group(3), match.group(4)
    new_record_id = basic_record_id.replace(
        '_c{}_t{}'.format(c, t),
        '_t{}'.format(rate_str)
    )

    return '{}.fasta'.format(new_record_id)
# end def


def group_repeat_records(repeat_records):
    grouped_repeat_records = dict()
    for rec in repeat_records:
        match = RECORD_ID_RE.match(str(rec.id))
        if match is None:
            sys.stderr.write('ERROR: rec.id `{}` is not not of expected format')
            sys.exit(1)
        # end if
        repeat_len = match.group(1)
        replicate_idx = match.group(2)
        group_key = (repeat_len, replicate_idx,)
        if group_key not in grouped_repeat_records:
            grouped_repeat_records[group_key] = {}
        # end if
        copy_idx = int(match.group(3))
        grouped_repeat_records[group_key][copy_idx] = rec
    # end for
    return grouped_repeat_records
# end def


def write_true_coords(true_repeat_coords, strand, true_repeats_bed_fpath):
    sep ='\t'
    with open(true_repeats_bed_fpath, 'wt') as out_handle:
        for chrom, start, end in true_repeat_coords:
            out_handle.write('{}\n'.format(sep.join([
                chrom,
                str(start),
                str(end),
                '.',
                strand,
            ])))
        # end for
    # end with
# end def

def parse_repat_len_from_repeat_id(repeat_id):
    return repeat_id.split('_')[2]
# end def


def insert_mock_repeat(chr_record,
                       seq_records_to_insert,
                       target_chrom_id='replicon_1'):
    # Use first record as representative for output naming
    repeat_descr = seq_records_to_insert[0].description

    orig_start_closed, orig_end_closed = get_orig_closed_coords(repeat_descr)

    seqs = [str(rec.seq) for rec in seq_records_to_insert]
    N = len(seqs)
    orig_chr_len = len(chr_record)

    curr_shoulder_len = SHOULDER_LEN
    if N > 1:
        max_repeat_len = max(map(len, seqs))
        curr_shoulder_len = SHOULDER_LEN + max_repeat_len
    # end if

    # Generate N + 1 positions (N insertion + original) with spacing ≥ SHOULDER_LEN
    coordinates_ok = False
    while not coordinates_ok:
        coords = [random.randint(0, orig_chr_len) for _ in range(N)]
        coords.sort()

        coordinates_ok = all(map(
            lambda coord: coord <= (orig_start_closed - curr_shoulder_len) \
                      or  coord >= (orig_end_closed   + curr_shoulder_len),
            coords
        ))
        if N > 1:
            coordinates_ok = coordinates_ok and all(map(
                lambda i: coords[i + 1] - coords[i] >= curr_shoulder_len,
                range(N)
            ))
        # end if
    # end while

    chr_str = str(chr_record.seq)
    orig_end_open = orig_end_closed + 1
    chr_str = mask_region(chr_str, orig_start_closed, orig_end_open)

    inserted_coords = []
    repeat_len_cumsum = 0
    orig_shift = 0

    for i, pos in enumerate(coords):
        seq = seqs[i]
        repeat_len = len(seq)
        adjusted_pos = pos + repeat_len_cumsum

        if pos < orig_start_closed:
            orig_shift += repeat_len
        # end if

        chr_str = chr_str[:adjusted_pos] + seq.lower() + chr_str[adjusted_pos:]
        inserted_coords.append((adjusted_pos, adjusted_pos + repeat_len,))
        repeat_len_cumsum += repeat_len
    # end for

    assert len(chr_str) == orig_chr_len + repeat_len_cumsum

    new_record = SeqRecord(Seq(chr_str), id=target_chrom_id, description='')

    orig_start_closed = orig_start_closed + orig_shift
    orig_end_open     = orig_end_closed + orig_shift + 1

    true_coords = [('replicon_1', orig_start_closed, orig_end_open)]
    for start, end in inserted_coords:
        true_coords.append((target_chrom_id, start, end))
    # end for

    strand = '+'
    if 'strand=-' in repeat_descr:
        strand = '-'
    # end if

    return new_record, true_coords, strand
# end def

def get_orig_closed_coords(orig_seq_description: str) -> tuple[int, int]:
    re_result = re.search(REPEAT_COORDS_RE, orig_seq_description)
    orig_start_closed = int(re_result.group(1))
    orig_end_closed   = int(re_result.group(2))
    return orig_start_closed, orig_end_closed
# end def

def mask_region(chr_str: str, start_closed: int, end_open: int) -> str:
    chr_len_before = len(chr_str)
    region_lower_seq = chr_str[start_closed : end_open].lower()
    chr_str = chr_str[:start_closed] + region_lower_seq + chr_str[end_open:]
    assert len(chr_str) == chr_len_before
    return chr_str
# end def

# >>> Proceed >>>

chr_record = next(iter(tuple(
    SeqIO.parse(IN_GENOME_FPATH, 'fasta')
)))
chr_record.id = 'replicon_1'

plasmid_record = next(iter(tuple(
    SeqIO.parse(IN_PLASMID_FPATH, 'fasta')
)))

chr_with_inserts_dir_path = os.path.join(
    IN_MOCK_REPEATS_DIRPATH,
    'chr_with_inserts'
)
if not os.path.isdir(chr_with_inserts_dir_path):
    os.makedirs(chr_with_inserts_dir_path)
# end if

replicate_id_list_fpath = os.path.join(
    IN_MOCK_REPEATS_DIRPATH,
    'replicate_id_list.txt'
)

true_repeat_coords_dirpath = os.path.join(
    IN_MOCK_REPEATS_DIRPATH,
    'true_repeat_locations',
)
if not os.path.isdir(true_repeat_coords_dirpath):
    os.makedirs(true_repeat_coords_dirpath)
# end if


sep = '\t'
rates = [None] + np.arange(RATE_FROM, RATE_TO+RATE_STEP, RATE_STEP).tolist()

with open(replicate_id_list_fpath, 'wt') as list_handle:

    for rate in rates:
        print('Doing repeats with rate = {}'.format(rate))
        repeat_fasta_fpath = make_mock_repat_fasta_fpath(rate)

        repeat_records = tuple(
            SeqIO.parse(repeat_fasta_fpath, 'fasta')
        )

        grouped_repeat_records = group_repeat_records(repeat_records)

        for group_key in sorted(grouped_repeat_records.keys()):
            dict_of_copies = grouped_repeat_records[group_key]
            sorted_copy_idxs = sorted(dict_of_copies.keys())

            # First copy → replicon_2 (plasmid)
            records_to_insert = [dict_of_copies[i] for i in sorted_copy_idxs[:1]]
            assert(len(records_to_insert)) == 1
            replicon_2_record_with_insert, true_coords_2, strand = insert_mock_repeat(
                plasmid_record,
                records_to_insert,
                'replicon_2'
            )
            # Remove original coords: they might be modified
            #   by the later insert_mock_repeat() call
            true_coords_2 = true_coords_2[1:]

            # Remaining dict_of_copies → replicon_1 (genome)
            if N_REPEAT_COPIES_TO_INSERT > 1:
                # strand is the same as at the first insert_mock_repeat() call
                records_to_insert = [dict_of_copies[i] for i in sorted_copy_idxs[1:]]
                assert(len(records_to_insert)) == N_REPEAT_COPIES_TO_INSERT - 1
                replicon_1_record_with_insert, true_coords_1, _ = insert_mock_repeat(
                    chr_record,
                    records_to_insert,
                    'replicon_1'
                )
            else:
                replicon_1_record_with_insert = chr_record
                repeat_descr = dict_of_copies[sorted_copy_idxs[0]].description
                orig_start_closed, orig_end_closed = get_orig_closed_coords(repeat_descr)
                orig_end_open = orig_end_closed + 1
                true_coords_1 = [('replicon_1', orig_start_closed, orig_end_open,)]
            # end if

            # Use first record as representative for output naming
            repeat_record = dict_of_copies[sorted_copy_idxs[0]]

            chr_with_insert_fpath = os.path.join(
                chr_with_inserts_dir_path,
                make_chr_with_insert_base_name(repeat_record, rate)
            )

            with open(chr_with_insert_fpath, 'wt') as tmp_out_handle:
                SeqIO.write(
                    [replicon_1_record_with_insert, replicon_2_record_with_insert,],
                    tmp_out_handle,
                    'fasta'
                )
            # end with

            curr_replicate_id = os.path.basename(chr_with_insert_fpath).replace(
                '.fasta',
                ''
            )
            list_handle.write('{}\n'.format(curr_replicate_id))

            true_coords_all = true_coords_1 + true_coords_2
            assert len(true_coords_all) == N_DESIRED_REPEAT_COPIES, \
                '{} != {}'.format(len(true_coords_all), N_DESIRED_REPEAT_COPIES)

            true_coord_bed_fpath = os.path.join(
                true_repeat_coords_dirpath,
                '{}.bed'.format(curr_replicate_id)
            )
            write_true_coords(true_coords_all, strand, true_coord_bed_fpath)
        # end for
    # end for
# end with


sys.exit(0)
