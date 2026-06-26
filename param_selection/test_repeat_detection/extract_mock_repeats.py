#!/usr/bin/env python3

import os
import sys
import random

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from mock_repeats_settings import MOCK_REPEAT_LENGTHS, \
                                  N_REPEAT_REPLICATES, \
                                  RANDOM_SEED


random.seed(RANDOM_SEED)

infpath = os.path.abspath(sys.argv[1])
outfpath = os.path.abspath(sys.argv[2])


chr_record = next(iter(tuple(
    SeqIO.parse(infpath, 'fasta')
)))
chr_len = len(chr_record)
chr_str = str(chr_record.seq)
copy_idx = '0'
rate = '0.00'

with open(outfpath, 'wt') as out_handle:
    for repeat_len in MOCK_REPEAT_LENGTHS:
        for replicate_idx in range(N_REPEAT_REPLICATES):

            start_coord = random.randint(0, (chr_len//2)-1)
            end_coord_open = start_coord + repeat_len
            strand = random.randint(0, 1)

            seq = Seq(
                chr_str[start_coord : end_coord_open]
            )
            if strand == 1:
                seq = seq.reverse_complement()
            # end if

            out_record = SeqRecord(
                seq,
                id='MR_l{}_r{}_c{}_t{}'.format(repeat_len, replicate_idx, copy_idx, rate),
                description='start={} end={} strand={}; all zero-based, closed'.format(
                    start_coord,
                    end_coord_open - 1,
                    '+' if strand == 0 else '-'
                )
            )

            SeqIO.write([out_record,], out_handle, 'fasta')
        # end for
    # end for
# end with

sys.exit(0)
