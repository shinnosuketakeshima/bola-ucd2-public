from scripts.codex_reanalysis.compare_dq_junctions import (
    difference_bins,
    junction_structure,
    overlap_loci,
    parse_hits,
    reverse_complement,
)

LEFT, RIGHT = 25581631, 25658061


def hit(q1, q2, r1, r2, identity=98.0, strand="+"):
    return {"q_start": q1, "q_end": q2, "r_start": r1, "r_end": r2, "strand": strand, "identity": identity,
            "length": q2 - q1 + 1}


def test_hits_are_shifted_and_oriented():
    (h,) = parse_hits("10\t20\t300\t290\t97.5\t11\n", 25569999)
    assert (h["r_start"], h["r_end"], h["strand"]) == (25570289, 25570299, "-")


def test_boundary_overlap_measures_breakpoint_uncertainty():
    # flank 3000; left anchor ends at query 3001, internal match starts 933 bp earlier (shared repeat),
    # internal match ends 15 bp after the right anchor starts (microhomology)
    hits = [
        hit(1, 3001, LEFT - 3000, LEFT),
        hit(2069, 6563, 25647043, 25651540, 96.0),
        hit(6562, 7451, 25651666, 25652553, 97.0),  # adjacent piece shorter than the seed minimum
        hit(4135, 5640, 25593194, 25594710, 92.0),  # repeat match elsewhere in the deleted interval
        hit(7437, 10437, RIGHT, RIGHT + 3000),
        hit(2069, 3000, 25640000, 25640931, 92.0),  # repeat copy shorter than the internal minimum
    ]
    row = junction_structure(hits, 10437, 3000, LEFT, RIGHT, 90.0, 1000)
    assert row["left_anchor_reference_end"] == LEFT and row["right_anchor_reference_start"] == RIGHT
    assert (row["internal_reference_start"], row["internal_reference_end"]) == (25647043, 25652553)
    assert row["internal_hits"] == 2
    assert row["left_boundary_overlap_bp"] == 933
    assert row["right_boundary_overlap_bp"] == 15
    assert overlap_loci(hits, row["left_overlap_query"]) == 3


def test_difference_bins_place_an_indel():
    source = "ACGT" * 250
    other = source[:600] + "T" * 20 + source[600:]
    bins = difference_bins(source, other, 500)
    assert bins[0] == 0 and bins[1] == 20


def test_reverse_complement():
    assert reverse_complement("AACGTN") == "NACGTT"
