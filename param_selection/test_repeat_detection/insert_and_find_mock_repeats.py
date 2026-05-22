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
SHOULDER_LEN = 100 # bp, around the original repeat


IN_GENOME_FPATH = os.path.abspath(sys.argv[1])
IN_MOCK_REPEATS_FILE = os.path.abspath(sys.argv[2])
IN_PIPELINE_DATA_DIR = os.path.abspath(sys.argv[3])
REVERLOR_FIND_FPATH = os.path.abspath(sys.argv[4])
TMP_DIRPATH = os.path.abspath(sys.argv[5])
OUT_FIND_REPEATS_ROOT_DIRPATH = os.path.abspath(sys.argv[6])



def make_mock_repat_fasta_fpath(rate=None):
    if rate is None:
        return IN_MOCK_REPEATS_FILE
    # end if

    rate_str = '{:.2f}'.format(rate).replace('.', 'dot')
    out_basename = os.path.basename(IN_MOCK_REPEATS_FILE).replace(
        '.fasta',
        '_{}_ms.fasta'.format(rate_str)
    )
    return os.path.join(
        IN_PIPELINE_DATA_DIR,
        out_basename
    )
# end def


def insert_mock_repeat(chr_record, mock_repeat_seq_record):
    chr_len = len(chr_record)
    chr_str = str(chr_record.seq)

    repeat_descr = mock_repeat_seq_record.description
    re_result = re.search(REPEAT_COORDS_RE, repeat_descr)

    start_coord = int(re_result.group(1))
    end_coord = int(re_result.group(2))

    start_coord_broader = start_coord - SHOULDER_LEN
    end_coord_broader = end_coord + SHOULDER_LEN

    insert_start_coord = random.randint(0, chr_len)
    # TODO: remove
    # while insert_start_coord >= start_coord_broader and insert_start_coord <= end_coord_broader:
    while insert_start_coord <= end_coord_broader:
        insert_start_coord = random.randint(0, chr_len)
    # end while

    end_coord_orig_open = end_coord + 1
    end_coord_inserted_open = insert_start_coord + len(mock_repeat_seq_record)

    true_repeat_coords = (
        (chr_record.id, start_coord, end_coord_orig_open),
        (chr_record.id, insert_start_coord, end_coord_inserted_open),
    )

    repeat_str = str(mock_repeat_seq_record.seq)
    new_chr_str = chr_str[:insert_start_coord] + repeat_str + chr_str[insert_start_coord:]
    new_chr_record = SeqRecord(
        Seq(new_chr_str),
        id=chr_record.id
    )

    return new_chr_record, true_repeat_coords
# end def

def run_reverlor_find(input_fasta_fpath, out_dir_path):
    cmd = ' '.join([
        'python3', REVERLOR_FIND_FPATH,
        input_fasta_fpath,
        out_dir_path,
    ])
    returncode = os.system(cmd)
    if returncode != 0:
        sys.stderr.write('Error running reverlor_find.py')
        sys.exit(1)
    # end if

    return os.path.join(out_dir_path, 'repeats_final.bed')
# end def

def make_curr_outdir_path(repeat_record, rate=None):
    if rate is None:
        base_name = 'default_{}'.format(repeat_record.id)
    else:
        base_name = 'default_{}_{:.2f}'.format(repeat_record.id, rate)
    # end def
    return os.path.join(
        OUT_FIND_REPEATS_ROOT_DIRPATH,
        base_name
    )
# end def

def write_true_coords(true_repeat_coords, true_repeats_bed_fpath):
    with open(true_repeats_bed_fpath, 'wt') as out_handle:
        for chrom, start, end in true_repeat_coords:
            out_handle.write(
                '{}\t{}\t{}\n'.format(
                    chrom,
                    start,
                    end
                )
            )
        # end for
    # end with
# end def



# >>> Proceed >>>

chr_record = next(iter(tuple(
    SeqIO.parse(IN_GENOME_FPATH, 'fasta')
)))

tmp_fasta_with_insert_fpath = os.path.join(
    TMP_DIRPATH,
    'chr_with_insert.fasta'
)

rates = [None] + np.arange(RATE_FROM, RATE_TO+RATE_STEP, RATE_STEP).tolist()

for rate in rates:
    repeat_fasta_fpath = make_mock_repat_fasta_fpath(rate)

    repeat_records = tuple(
        SeqIO.parse(repeat_fasta_fpath, 'fasta')
    )
    
    for repeat_record in repeat_records:
        chr_record_with_insert, true_repeat_coords = insert_mock_repeat(
            chr_record,
            repeat_record
        )

        with open(tmp_fasta_with_insert_fpath, 'wt') as tmp_out_handle:
            SeqIO.write(
                [chr_record_with_insert,],
                tmp_out_handle,
                'fasta'
            )
        # end with

        curr_out_dir_path = make_curr_outdir_path(repeat_record, rate)
        out_repeats_fpath = run_reverlor_find(
            tmp_fasta_with_insert_fpath,
            curr_out_dir_path
        )

        true_repeats_bed_fpath = os.path.join(
            curr_out_dir_path,
            os.path.basename(out_repeats_fpath).replace('.bed', '.true.bed')
        )
        write_true_coords(true_repeat_coords, true_repeats_bed_fpath)
    # end for
# end for


if os.path.isfile(tmp_fasta_with_insert_fpath):
    os.unlink(tmp_fasta_with_insert_fpath)
# end if


sys.exit(0)
