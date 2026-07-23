import pytest

from reverlor.src.FindArgs import FindArgs
from reverlor.src.bed_lib import (
    RepeatRegion,
    VerifyResult,
    read_bed_to_regions,
    verify_results_to_bed,
    merge_features,
)


# >>> read_bed_to_regions tests >>>

def test_single_region(tmp_path):
    bed = tmp_path / 'one.bed'
    bed.write_text('chr1\t100\t200\n')
    regions = read_bed_to_regions(str(bed))
    assert len(regions) == 1
    assert regions[0].ref_id == 'chr1'
    assert regions[0].start == 100
    assert regions[0].end == 200
# end def


def test_multiple_regions(tmp_path):
    bed = tmp_path / 'multi.bed'
    bed.write_text('chr1\t100\t200\nchr2\t300\t400\nchr1\t500\t600\n')
    regions = read_bed_to_regions(str(bed))
    assert len(regions) == 3
    assert regions[0].ref_id == 'chr1'
    assert regions[0].start == 100
    assert regions[1].ref_id == 'chr2'
    assert regions[1].start == 300
    assert regions[2].ref_id == 'chr1'
    assert regions[2].start == 500
# end def


def test_empty_file(tmp_path):
    bed = tmp_path / 'empty.bed'
    bed.write_text('')
    regions = read_bed_to_regions(str(bed))
    assert regions == []
# end def


def test_whitespace_stripped(tmp_path):
    bed = tmp_path / 'ws.bed'
    bed.write_text('  chr1\t100  \t  200\nchr2 \t300  \t 400\n')
    regions = read_bed_to_regions(str(bed))
    assert len(regions) == 2
    assert regions[0].ref_id == 'chr1'
    assert regions[1].ref_id == 'chr2'
# end def


def test_coordinates_are_integers(tmp_path):
    bed = tmp_path / 'ints.bed'
    bed.write_text('chr1\t100\t200\n')
    regions = read_bed_to_regions(str(bed))
    assert isinstance(regions[0].start, int)
    assert isinstance(regions[0].end, int)
# end def


# <<< read_bed_to_regions tests <<<


# >>> verify_results_to_bed tests >>>

def test_single_result(tmp_path):
    out = tmp_path / 'out.bed'
    vr = VerifyResult(region=RepeatRegion('chr1', 100, 200), num_read_throughs=7)
    verify_results_to_bed([vr], str(out))
    lines = out.read_text().strip().split('\n')
    assert len(lines) == 1
    assert lines[0] == 'chr1\t100\t200\trepeat\t7'
# end def


def test_multiple_results(tmp_path):
    out = tmp_path / 'out.bed'
    results = [
        VerifyResult(region=RepeatRegion('chr1', 100, 200), num_read_throughs=7),
        VerifyResult(region=RepeatRegion('chr2', 300, 400), num_read_throughs=2),
        VerifyResult(region=RepeatRegion('chr1', 500, 600), num_read_throughs=15),
    ]
    verify_results_to_bed(results, str(out))
    lines = out.read_text().strip().split('\n')
    assert len(lines) == 3
    assert lines[0] == 'chr1\t100\t200\trepeat\t7'
    assert lines[1] == 'chr2\t300\t400\trepeat\t2'
    assert lines[2] == 'chr1\t500\t600\trepeat\t15'
# end def


def test_empty_list(tmp_path):
    out = tmp_path / 'out.bed'
    verify_results_to_bed([], str(out))
    assert out.exists()
    assert out.read_text() == ''
# end def


def test_output_is_tab_separated(tmp_path):
    out = tmp_path / 'out.bed'
    vr = VerifyResult(region=RepeatRegion('chr1', 100, 200), num_read_throughs=7)
    verify_results_to_bed([vr], str(out))
    line = out.read_text().strip()
    parts = line.split('\t')
    assert len(parts) == 5
    assert parts[0] == 'chr1'
    assert parts[3] == 'repeat'
# end def


def test_round_trip(tmp_path):
    bed = tmp_path / 'round_trip.bed'
    results = [
        VerifyResult(region=RepeatRegion('chrX', 10, 50), num_read_throughs=3),
        VerifyResult(region=RepeatRegion('chrY', 200, 800), num_read_throughs=12),
    ]
    verify_results_to_bed(results, str(bed))
    regions = read_bed_to_regions(str(bed))
    assert len(regions) == 2
    assert regions[0].ref_id == 'chrX'
    assert regions[0].start == 10
    assert regions[0].end == 50
    assert regions[1].ref_id == 'chrY'
    assert regions[1].start == 200
    assert regions[1].end == 800
# end def


# <<< verify_results_to_bed tests <<<


# >>> merge_features tests >>>

def _make_find_args(min_repeat_len, min_repeat_interval):
    return FindArgs(
        fasta_fpath='/dev/null',
        output_dir='/tmp',
        min_repeat_len=min_repeat_len,
        min_repeat_interval=min_repeat_interval,
    )
# end def


def _read_merged_bed(fpath):
    regions = []
    with open(fpath, 'rt') as fh:
        for line in fh:
            vals = line.strip().split('\t')
            if vals == ['']:
                continue
            # end if
            regions.append((vals[0], int(vals[1]), int(vals[2])))
        # end for
    # end with
    return regions
# end def


def test_single_region_passes_through(tmp_path):
    in_bed = tmp_path / 'in.bed'
    out_bed = tmp_path / 'out.bed'
    in_bed.write_text('chr1\t100\t400\n')
    args = _make_find_args(min_repeat_len=200, min_repeat_interval=100)
    merge_features(args, str(in_bed), str(out_bed))
    regions = _read_merged_bed(out_bed)
    assert len(regions) == 1
    assert regions[0] == ('chr1', 100, 400)
# end def


def test_overlapping_regions_merge(tmp_path):
    in_bed = tmp_path / 'in.bed'
    out_bed = tmp_path / 'out.bed'
    in_bed.write_text('chr1\t100\t300\nchr1\t200\t400\n')
    args = _make_find_args(min_repeat_len=50, min_repeat_interval=100)
    merge_features(args, str(in_bed), str(out_bed))
    regions = _read_merged_bed(out_bed)
    assert len(regions) == 1
    assert regions[0] == ('chr1', 100, 400)
# end def


def test_adjacent_within_interval_merge(tmp_path):
    in_bed = tmp_path / 'in.bed'
    out_bed = tmp_path / 'out.bed'
    in_bed.write_text('chr1\t100\t200\nchr1\t210\t300\n')
    args = _make_find_args(min_repeat_len=50, min_repeat_interval=20)
    merge_features(args, str(in_bed), str(out_bed))
    regions = _read_merged_bed(out_bed)
    assert len(regions) == 1
    assert regions[0] == ('chr1', 100, 300)
# end def


def test_adjacent_beyond_interval_stay_separate(tmp_path):
    in_bed = tmp_path / 'in.bed'
    out_bed = tmp_path / 'out.bed'
    in_bed.write_text('chr1\t100\t200\nchr1\t350\t450\n')
    args = _make_find_args(min_repeat_len=50, min_repeat_interval=100)
    merge_features(args, str(in_bed), str(out_bed))
    regions = _read_merged_bed(out_bed)
    assert len(regions) == 2
    assert regions[0] == ('chr1', 100, 200)
    assert regions[1] == ('chr1', 350, 450)
# end def


def test_short_regions_filtered_out(tmp_path):
    in_bed = tmp_path / 'in.bed'
    out_bed = tmp_path / 'out.bed'
    in_bed.write_text('chr1\t100\t200\nchr1\t500\t900\n')
    args = _make_find_args(min_repeat_len=200, min_repeat_interval=100)
    merge_features(args, str(in_bed), str(out_bed))
    regions = _read_merged_bed(out_bed)
    assert len(regions) == 1
    assert regions[0] == ('chr1', 500, 900)
# end def


def test_empty_file(tmp_path):
    in_bed = tmp_path / 'in.bed'
    out_bed = tmp_path / 'out.bed'
    in_bed.write_text('')
    args = _make_find_args(min_repeat_len=200, min_repeat_interval=100)
    merge_features(args, str(in_bed), str(out_bed))
    regions = _read_merged_bed(out_bed)
    assert regions == []
# end def


def test_unsorted_input_gets_sorted(tmp_path):
    in_bed = tmp_path / 'in.bed'
    out_bed = tmp_path / 'out.bed'
    in_bed.write_text('chr1\t500\t600\nchr1\t100\t200\n')
    args = _make_find_args(min_repeat_len=50, min_repeat_interval=100)
    merge_features(args, str(in_bed), str(out_bed))
    regions = _read_merged_bed(out_bed)
    assert len(regions) == 2
    assert regions[0] == ('chr1', 100, 200)
    assert regions[1] == ('chr1', 500, 600)
# end def


def test_merge_then_filter(tmp_path):
    in_bed = tmp_path / 'in.bed'
    out_bed = tmp_path / 'out.bed'
    in_bed.write_text('chr1\t100\t150\nchr1\t160\t200\n')
    args = _make_find_args(min_repeat_len=100, min_repeat_interval=20)
    merge_features(args, str(in_bed), str(out_bed))
    regions = _read_merged_bed(out_bed)
    assert len(regions) == 1
    assert regions[0] == ('chr1', 100, 200)
# end def


# <<< merge_features tests <<<
