import os
import tempfile
from pathlib import Path

import pytest
from Bio import SeqIO

from reverlor.src.FindArgs import FindArgs
from reverlor.src.bed_lib import read_bed_to_regions
from reverlor.src.find_repeats_minimap2 import find_repeats


DATA_DIR = Path(__file__).parent / 'data' / 'for_find_repeats'
REPEAT_LEN = 256
MAX_ERROR_BP = 5
INSERTS_CHR1 = [1000, 3000, 5000]
INSERTS_CHR2 = [2000]


# >>> Helper functions >>>

def _load_refs(paths):
    refs = {}
    for path in paths:
        for rec in SeqIO.parse(str(path), 'fasta'):
            refs[rec.id] = str(rec.seq)
        # end for
    # end for
    return refs
# end def


def _insert_repeats(seqs, insertions, repeat_seq):
    coords_by_ref = {}
    for ref_id, pos in insertions:
        coords_by_ref.setdefault(ref_id, []).append(pos)
    # end for

    result = dict(seqs)
    for ref_id, positions in coords_by_ref.items():
        seq = result[ref_id]
        for pos in sorted(positions, reverse=True):
            seq = seq[:pos] + repeat_seq + seq[pos:]
        # end for
        result[ref_id] = seq
    # end for
    return result
# end def


def _write_fasta(seqs, out_path):
    with open(out_path, 'w') as fh:
        for name, seq in seqs.items():
            fh.write(f'>{name}\n')
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + '\n')
            # end for
        # end for
    # end with
# end def


def _compute_true_positions(insertions, repeat_len):
    coords_by_ref = {}
    for ref_id, pos in insertions:
        coords_by_ref.setdefault(ref_id, []).append(pos)
    # end for

    positions = []
    for ref_id, coords in coords_by_ref.items():
        for i, pos in enumerate(sorted(coords)):
            start = pos + i * repeat_len
            end = start + repeat_len
            positions.append((ref_id, start, end))
        # end for
    # end for
    return positions
# end def


def _assert_results(bed_path, insertions, repeat_len, expected_count):
    regions = read_bed_to_regions(bed_path)
    assert len(regions) == expected_count, (
        f'Expected {expected_count} regions, got {len(regions)}: '
        f'{[(r.ref_id, r.start, r.end) for r in regions]}'
    )
    # end if

    if expected_count == 0:
        return
    # end if

    true_positions = _compute_true_positions(insertions, repeat_len)

    true_coords_by_ref = {}
    for ref_id, start, end in true_positions:
        true_coords_by_ref.setdefault(ref_id, []).append((start, end))
    # end for

    detected_coords_by_ref = {}
    for r in regions:
        detected_coords_by_ref.setdefault(r.ref_id, []).append((r.start, r.end))
    # end for

    assert set(detected_coords_by_ref.keys()) == set(true_coords_by_ref.keys()), (
        f'Chromosome mismatch: detected {set(detected_coords_by_ref.keys())}, '
        f'expected {set(true_coords_by_ref.keys())}'
    )
    # end if

    for ref_id in true_coords_by_ref:
        det = sorted(detected_coords_by_ref[ref_id])
        tru = sorted(true_coords_by_ref[ref_id])
        assert len(det) == len(tru), (
            f'{ref_id}: expected {len(tru)} regions, got {len(det)}'
        )
        # end if

        for (d_start, d_end), (t_start, t_end) in zip(det, tru):
            assert abs(d_start - t_start) <= MAX_ERROR_BP, (
                f'{ref_id}: detected start {d_start} too far from true {t_start} '
                f'(error {abs(d_start - t_start)} > {MAX_ERROR_BP})'
            )
            # end if
            assert abs(d_end - t_end) <= MAX_ERROR_BP, (
                f'{ref_id}: detected end {d_end} too far from true {t_end} '
                f'(error {abs(d_end - t_end)} > {MAX_ERROR_BP})'
            )
            # end if
        # end for
    # end for
# end def


def _run_find_repeats(ref_paths, insertions, min_repeat_len=200):
    repeat_seq = str(list(SeqIO.parse(
        str(DATA_DIR / 'repeat.fasta'), 'fasta'
    ))[0].seq)

    seqs = _load_refs(ref_paths)
    seqs = _insert_repeats(seqs, insertions, repeat_seq)

    fasta_path = tempfile.mktemp(suffix='.fasta')
    _write_fasta(seqs, fasta_path)

    out_dir = tempfile.mkdtemp()
    args = FindArgs(
        fasta_fpath=fasta_path,
        output_dir=out_dir,
        min_repeat_len=min_repeat_len,
    )
    bed_path = find_repeats(args)
    os.unlink(fasta_path)
    return bed_path
# end def


# <<< Helper functions <<<


# >>> Core tests >>>

def test_no_inserts(tmp_path):
    bed_path = _run_find_repeats(
        ref_paths=[
            DATA_DIR / 'some_seq_no_repeats_1.fasta',
            DATA_DIR / 'some_seq_no_repeats_2.fasta',
        ],
        insertions=[],
    )
    _assert_results(bed_path, insertions=[], repeat_len=REPEAT_LEN, expected_count=0)
# end def


def test_single_copy_no_partner():
    bed_path = _run_find_repeats(
        ref_paths=[DATA_DIR / 'some_seq_no_repeats_1.fasta'],
        insertions=[('chr1', 1000)],
    )
    _assert_results(
        bed_path,
        insertions=[('chr1', 1000)],
        repeat_len=REPEAT_LEN,
        expected_count=0,
    )
# end def


def test_two_copies_one_ref():
    insertions = [
        ('chr1', INSERTS_CHR1[0]),
        ('chr1', INSERTS_CHR1[1]),
    ]
    bed_path = _run_find_repeats(
        ref_paths=[DATA_DIR / 'some_seq_no_repeats_1.fasta'],
        insertions=insertions,
    )
    _assert_results(bed_path, insertions=insertions, repeat_len=REPEAT_LEN, expected_count=2)
# end def


def test_two_copies_min_len_300():
    insertions = [
        ('chr1', INSERTS_CHR1[0]),
        ('chr1', INSERTS_CHR1[1]),
    ]
    bed_path = _run_find_repeats(
        ref_paths=[DATA_DIR / 'some_seq_no_repeats_1.fasta'],
        insertions=insertions,
        min_repeat_len=300,
    )
    _assert_results(bed_path, insertions=insertions, repeat_len=REPEAT_LEN, expected_count=0)
# end def


def test_three_copies_one_ref():
    insertions = [
        ('chr1', INSERTS_CHR1[0]),
        ('chr1', INSERTS_CHR1[1]),
        ('chr1', INSERTS_CHR1[2]),
    ]
    bed_path = _run_find_repeats(
        ref_paths=[DATA_DIR / 'some_seq_no_repeats_1.fasta'],
        insertions=insertions,
    )
    _assert_results(bed_path, insertions=insertions, repeat_len=REPEAT_LEN, expected_count=3)
# end def


def test_one_copy_per_ref():
    insertions = [
        ('chr1', INSERTS_CHR1[0]),
        ('chr2', INSERTS_CHR2[0]),
    ]
    bed_path = _run_find_repeats(
        ref_paths=[
            DATA_DIR / 'some_seq_no_repeats_1.fasta',
            DATA_DIR / 'some_seq_no_repeats_2.fasta',
        ],
        insertions=insertions,
    )
    _assert_results(bed_path, insertions=insertions, repeat_len=REPEAT_LEN, expected_count=2)
# end def


def test_two_copies_chr1_one_copy_chr2():
    insertions = [
        ('chr1', INSERTS_CHR1[0]),
        ('chr1', INSERTS_CHR1[1]),
        ('chr2', INSERTS_CHR2[0]),
    ]
    bed_path = _run_find_repeats(
        ref_paths=[
            DATA_DIR / 'some_seq_no_repeats_1.fasta',
            DATA_DIR / 'some_seq_no_repeats_2.fasta',
        ],
        insertions=insertions,
    )
    _assert_results(bed_path, insertions=insertions, repeat_len=REPEAT_LEN, expected_count=3)
# end def


# <<< Core tests <<<


# >>> Edge-case tests >>>

def test_repeat_at_start_of_sequence():
    insertions = [('chr1', 0), ('chr1', 3000)]
    bed_path = _run_find_repeats(
        ref_paths=[DATA_DIR / 'some_seq_no_repeats_1.fasta'],
        insertions=insertions,
    )
    _assert_results(bed_path, insertions=insertions, repeat_len=REPEAT_LEN, expected_count=2)
# end def


def test_repeat_near_end_of_sequence():
    insertions = [('chr1', 3000), ('chr1', 6700)]
    bed_path = _run_find_repeats(
        ref_paths=[DATA_DIR / 'some_seq_no_repeats_1.fasta'],
        insertions=insertions,
    )
    _assert_results(bed_path, insertions=insertions, repeat_len=REPEAT_LEN, expected_count=2)
# end def


# <<< Edge-case tests <<<
