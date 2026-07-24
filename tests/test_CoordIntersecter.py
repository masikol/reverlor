import pathlib

import pytest

from reverlor.src.CoordIntersecter import CoordIntersecter


BAM_PATH = str(
    pathlib.Path(__file__).parent / 'data' / 'for_CoordIntersecter' / 'some_mapping.srt.bam'
)


@pytest.fixture
def ci():
    return CoordIntersecter(BAM_PATH)
# end def


# >>> intersect_coords tests >>>

def test_both_positions_beyond_coverage(ci):
    assert ci.intersect_coords('Rho_K11-5_v1_chr', 1000, 1200) == frozenset()
# end def


def test_one_position_beyond_coverage(ci):
    assert ci.intersect_coords('Rho_K11-5_v1_chr', 50, 1200) == frozenset()
# end def


def test_one_read_spans_both(ci):
    result = ci.intersect_coords('Rho_K11-5_v1_chr', 660, 860)
    assert result == frozenset({'Rho_K11-5_v1_chr_2'})
# end def


def test_multiple_cover_pos1_one_covers_pos2(ci):
    result = ci.intersect_coords('Rho_K11-5_v1_chr', 489, 800)
    assert result == frozenset({'Rho_K11-5_v1_chr_2'})
# end def


def test_three_reads_span_both(ci):
    result = ci.intersect_coords('Rho_K11-5_v1_chr', 245, 445)
    assert result == frozenset({
        'Rho_K11-5_v1_chr_3',
        'Rho_K11-5_v1_chr_4',
        'Rho_K11-5_v1_chr_5',
    })
# end def


def test_intersection_filters_correctly(ci):
    result = ci.intersect_coords('Rho_K11-5_v1_chr', 50, 600)
    assert result == frozenset({'Rho_K11-5_v1_chr_4'})
# end def


def test_single_position_returns_all_covering_reads(ci):
    result = ci.intersect_coords('Rho_K11-5_v1_chr', 605, 605)
    assert result == frozenset({
        'Rho_K11-5_v1_chr_2',
        'Rho_K11-5_v1_chr_4',
    })
# end def


# <<< intersect_coords tests <<<


# >>> edge cases >>>

def test_negative_pos_1_raises(ci):
    with pytest.raises(ValueError, match='pos_1'):
        ci.intersect_coords('Rho_K11-5_v1_chr', -1, 100)
    # end with
# end def


def test_negative_pos_2_raises(ci):
    with pytest.raises(ValueError, match='pos_2'):
        ci.intersect_coords('Rho_K11-5_v1_chr', 100, -1)
    # end with
# end def


def test_nonexistent_rname_raises(ci):
    with pytest.raises(ValueError):
        ci.intersect_coords('nonexistent', 0, 100)
    # end with
# end def


# <<< edge cases <<<
