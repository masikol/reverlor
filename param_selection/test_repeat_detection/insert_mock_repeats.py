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
SHOULDER_LEN = 200 # bp, around the original repeat

IN_GENOME_FPATH = os.path.abspath(sys.argv[1])
IN_MOCK_REPEATS_DIRPATH = os.path.abspath(sys.argv[2])
MASK_ALL_BUT_REPEATS: bool = sys.argv[3] == '1'

IN_MOCK_REPEATS_FPATH = os.path.join(
    IN_MOCK_REPEATS_DIRPATH,
    'mock_repeats.fasta'
)



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
    record_id = str(mock_repeat_seq_record.id)
    if rate is None:
        return '{}.fasta'.format(record_id)
    # end if
    return '{}_{:.2f}.fasta'.format(record_id, rate)
# end def


def insert_mock_repeat(chr_record,
                       mock_repeat_seq_record,
                       mask_all_but_repeats=False):
    chr_len = len(chr_record)
    chr_str = str(chr_record.seq)

    repeat_descr = mock_repeat_seq_record.description
    re_result = re.search(REPEAT_COORDS_RE, repeat_descr)

    orig_start_coord_closed = int(re_result.group(1))
    orig_end_coord_closed = int(re_result.group(2))

    end_coord_broader = orig_end_coord_closed + SHOULDER_LEN
    insert_start_coord_closed = random.randint(end_coord_broader, chr_len)
    insert_end_coord_open = insert_start_coord_closed + len(mock_repeat_seq_record)

    orig_end_coord_open = orig_end_coord_closed + 1
    true_repeat_coords = (
        (chr_record.id, orig_start_coord_closed, orig_end_coord_open),
        (chr_record.id, insert_start_coord_closed, insert_end_coord_open),
    )

    repeat_str = str(mock_repeat_seq_record.seq)
    end_coord_inserted_closed = insert_end_coord_open - 1
    new_chr_str = chr_str[:insert_start_coord_closed] + repeat_str + chr_str[insert_start_coord_closed:]

    if mask_all_but_repeats:
        new_chr_str = mask_all_but_repeats(
            new_chr_str,
            orig_start_coord_closed,
            orig_end_coord_closed,
            insert_start_coord_closed,
            end_coord_inserted_closed
        )
    # end if

    new_chr_record = SeqRecord(
        Seq(new_chr_str),
        id=chr_record.id,
        description=''
    )

    strand = '+'
    if 'strand=-' in repeat_descr:
        strand = '-'
    # end if

    return new_chr_record, true_repeat_coords, strand
# end def

def mask_all_but_repeats(chr_str,
                         orig_start_coord_closed,
                         orig_end_coord_closed,
                         insert_start_coord_closed,
                         end_coord_inserted_closed):
    new_chr_chars = list(chr_str)

    for i in range(0, orig_start_coord_closed):
        new_chr_chars[i] = 'n'
    # end for
    for i in range(orig_end_coord_closed+1, insert_start_coord_closed):
        new_chr_chars[i] = 'n'
    # end for
    for i in range(end_coord_inserted_closed+1, len(new_chr_chars)):
        new_chr_chars[i] = 'n'
    # end for

    return ''.join(new_chr_chars)
# end def

def make_curr_outdir_path(repeat_record, rate=None):
    if rate is None:
        base_name = 'default_{}'.format(repeat_record.id)
    else:
        base_name = 'default_{}_{:.2f}'.format(repeat_record.id, rate)
    # end def
    return os.path.join(
        IN_MOCK_REPEATS_DIRPATH,
        'true_repeat_locations',
        base_name
    )
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



# >>> Proceed >>>

chr_record = next(iter(tuple(
    SeqIO.parse(IN_GENOME_FPATH, 'fasta')
)))

chr_with_inserts_dir_path = os.path.join(
    IN_MOCK_REPEATS_DIRPATH,
    'chr_with_inserts'
)
if not os.path.isdir(chr_with_inserts_dir_path):
    os.makedirs(chr_with_inserts_dir_path)
# end if

test_combintations_fpath = os.path.join(
    IN_MOCK_REPEATS_DIRPATH,
    'test_combintations.tsv'
)


sep = '\t'
rates = [None] + np.arange(RATE_FROM, RATE_TO+RATE_STEP, RATE_STEP).tolist()

with open(test_combintations_fpath, 'wt') as combin_handle:

    combin_handle.write('{}\n'.format(sep.join([
        'repeat_len',
        'repeat_id',
        'rate',
        'chr_with_insert_fname',
        'true_repeats_bed_fname',
        'pred_repeats_bed_fname',
    ])))

    for rate in rates:
        repeat_fasta_fpath = make_mock_repat_fasta_fpath(rate)

        repeat_records = tuple(
            SeqIO.parse(repeat_fasta_fpath, 'fasta')
        )
        
        for repeat_record in repeat_records:
            chr_record_with_insert, true_repeat_coords, strand = insert_mock_repeat(
                chr_record,
                repeat_record,
                MASK_ALL_BUT_REPEATS
            )

            chr_with_insert_fpath = os.path.join(
                chr_with_inserts_dir_path,
                make_chr_with_insert_base_name(repeat_record, rate)
            )

            with open(chr_with_insert_fpath, 'wt') as tmp_out_handle:
                SeqIO.write(
                    [chr_record_with_insert,],
                    tmp_out_handle,
                    'fasta'
                )
            # end with

            curr_out_dir_path = make_curr_outdir_path(repeat_record, rate)
            if not os.path.isdir(curr_out_dir_path):
                os.makedirs(curr_out_dir_path)
            # end if

            pred_repeats_bed_fpath = os.path.join(
                curr_out_dir_path,
                'repeats_final.bed'
            )

            true_repeats_bed_fpath = os.path.join(
                curr_out_dir_path,
                os.path.basename(pred_repeats_bed_fpath).replace('.bed', '.true.bed')
            )
            write_true_coords(true_repeat_coords, strand, true_repeats_bed_fpath)

            rate_str = '0.00'
            if rate is not None:
                rate_str = '{:.2f}'.format(rate)
            # end if

            combin_handle.write('{}\n'.format(sep.join([
                parse_repat_len_from_repeat_id(repeat_record.id),
                str(repeat_record.id),
                rate_str,
                os.path.basename(chr_with_insert_fpath),
                os.path.basename(true_repeats_bed_fpath),
                os.path.basename(pred_repeats_bed_fpath),
            ])))
        # end for
    # end for
# end with


sys.exit(0)
