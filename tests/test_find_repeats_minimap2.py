import os
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio import SeqIO

from reverlor.src.FindArgs import FindArgs
from reverlor.src.ReverlorArgs import ReverlorArgs
from reverlor.src.bed_lib import read_bed_to_regions
from reverlor.src.find_repeats_minimap2 import (
    _filter_hit,
    _filter_hit_by_pident,
    find_repeats,
)


DATA_DIR = Path(__file__).parent / 'data' / 'for_find_repeats'
REPEAT_LEN = 256
MAX_ERROR_BP = 5
INSERTS_CHR1 = [1000, 3000, 5000]
INSERTS_CHR2 = [2000]
FLOAT_EPSILON = 1e-6

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


# >>> min_pident tests >>>

def test_filter_hit_accepts_pident_at_threshold():
    args = FindArgs(
        fasta_fpath='/dev/null',
        output_dir='/tmp',
        min_pident=0.9,
    )
    hit = SimpleNamespace(mlen=90, blen=100)
    assert _filter_hit_by_pident(hit, args)
    assert _filter_hit(hit, args)
# end def


def test_filter_hit_rejects_pident_below_threshold():
    args = FindArgs(
        fasta_fpath='/dev/null',
        output_dir='/tmp',
        min_pident=0.9,
    )
    hit = SimpleNamespace(mlen=89, blen=100)
    assert not _filter_hit_by_pident(hit, args)
    assert not _filter_hit(hit, args)
# end def


def test_filter_hit_rejects_zero_blen():
    args = FindArgs(
        fasta_fpath='/dev/null',
        output_dir='/tmp',
        min_pident=0.0,
    )
    assert not _filter_hit_by_pident(SimpleNamespace(mlen=0, blen=0), args)
# end def


def test_min_pident_defaults_to_zero():
    find_args = FindArgs(fasta_fpath='/dev/null', output_dir='/tmp')
    reverlor_args = ReverlorArgs(
        fasta_fpath='/dev/null',
        input_bam_fpath='/dev/null',
        output_dir='/tmp',
    )
    assert abs(find_args.min_pident - 0.0) < FLOAT_EPSILON
    assert abs(reverlor_args.min_pident - 0.0) < FLOAT_EPSILON
# end def


def test_find_args_converts_min_pident_percent_to_ratio(monkeypatch, tmp_path):
    fasta_path = tmp_path / 'input.fasta'
    fasta_path.touch()
    monkeypatch.setattr(sys, 'argv', [
        'reverlor_find',
        str(fasta_path),
        str(tmp_path / 'out'),
        '--min-pident',
        '95.0',
    ])

    args = FindArgs.parse_args()

    assert abs(args.min_pident - 0.95) < FLOAT_EPSILON
# end def


def test_reverlor_args_converts_min_pident_percent_to_ratio(monkeypatch, tmp_path):
    fasta_path = tmp_path / 'input.fasta'
    fasta_path.touch()
    bam_path = tmp_path / 'input.bam'
    bam_path.touch()
    monkeypatch.setattr(sys, 'argv', [
        'reverlor',
        str(fasta_path),
        str(bam_path),
        str(tmp_path / 'out'),
        '--min-pident',
        '100.0',
    ])

    args = ReverlorArgs.parse_args()

    assert abs(args.min_pident - 1.0) < FLOAT_EPSILON
# end def


def test_find_args_propagates_min_pident_ratio():
    reverlor_args = ReverlorArgs(
        fasta_fpath='/dev/null',
        input_bam_fpath='/dev/null',
        output_dir='/tmp',
        min_pident=0.85,
    )

    find_args = FindArgs.from_reverlor_args(reverlor_args)

    assert abs(find_args.min_pident - 0.85) < FLOAT_EPSILON
# end def


def test_find_args_rejects_min_pident_above_100(monkeypatch, tmp_path):
    fasta_path = tmp_path / 'input.fasta'
    fasta_path.touch()
    monkeypatch.setattr(sys, 'argv', [
        'reverlor_find',
        str(fasta_path),
        str(tmp_path / 'out'),
        '--min-pident',
        '100.1',
    ])

    with pytest.raises(SystemExit) as exc_info:
        FindArgs.parse_args()
    # end with

    assert exc_info.value.code == 1
# end def


def test_reverlor_args_rejects_negative_min_pident(monkeypatch, tmp_path):
    fasta_path = tmp_path / 'input.fasta'
    fasta_path.touch()
    bam_path = tmp_path / 'input.bam'
    bam_path.touch()
    monkeypatch.setattr(sys, 'argv', [
        'reverlor',
        str(fasta_path),
        str(bam_path),
        str(tmp_path / 'out'),
        '--min-pident',
        '-0.1',
    ])

    with pytest.raises(SystemExit) as exc_info:
        ReverlorArgs.parse_args()
    # end with

    assert exc_info.value.code == 1
# end def


# <<< min_pident tests <<<


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
