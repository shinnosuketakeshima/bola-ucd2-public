from scripts.codex_reanalysis.compare_public_dq_configurations import (
    classify,
    diagnostic_windows,
    locate_span,
    paf_blocks,
)


def paf(t_start, t_end, strand, matches, block, q_start=0, q_end=None):
    q_end = q_end if q_end is not None else block
    return "\t".join(map(str, ["q", 200000, q_start, q_end, strand, "chr", 5e7, t_start, t_end, matches, block, 60]))


def test_locate_span_ignores_distant_paralog_hits():
    text = "\n".join([
        paf(26_000_000, 26_150_000, "+", 140000, 150000),
        paf(26_160_000, 26_190_000, "+", 28000, 30000),
        paf(31_000_000, 31_010_000, "-", 9000, 10000),  # >1 Mb away
        paf(26_200_000, 26_201_000, "+", 900, 1000),  # below the block minimum
    ])
    assert locate_span(text, 5000) == (26_000_001, 26_190_000, "+")
    assert locate_span("", 5000) is None


def test_paf_blocks_are_one_based_with_identity():
    (block,) = paf_blocks(paf(99, 1099, "-", 980, 1000, 9, 1009), 500)
    assert (block["target_start"], block["target_end"], block["query_start"]) == (100, 1099, 10)
    assert block["identity"] == 0.98 and block["strand"] == "-"
    assert paf_blocks(paf(99, 1099, "-", 980, 1000), 2000) == []


def test_diagnostic_windows_use_each_target_frame():
    status = {
        "selected_deletion": {"contig_left": 91774, "contig_right": 96310, "reference_left": 25581631},
        "selected_insertion": {"contig_left": 72164, "contig_right": 97927, "reference_left": 25576670,
                               "reference_right": 25576909},
        "parameters": {"insertion_source": "ID06:contig_9"},
    }
    hits = [{"sample_id": "ID06", "contig": "contig_9", "family": "DQA", "attributed_locus": "DQA2",
             "contig_start": "104063", "contig_end": "109598"}]
    features = {"BOLA-DQA2": (25583031, 25589201), "BOLA-DQB": (25607507, 25620857),
                "BOLA-DQA5": (25631958, 25643878)}
    windows = diagnostic_windows(status, hits, features, 25450000)
    assert windows["deletion_contig"] == [("deletion_junction", 91774, 96310)]
    assert ("insertion_config_DQA2", 104063, 109598) in windows["insertion_contig"]
    assert ("reference_deletion_point", 131631, 131631) in windows["reference"]
    assert ("reference_DQB_DQA5", 157507, 193878) in windows["reference"]


def test_classification_needs_every_diagnostic_of_a_configuration():
    assert classify({"deletion_junction": True}) == "A"
    assert classify({"insertion_left_edge": True, "insertion_whole_segment": True,
                     "insertion_right_edge": True}) == "B"
    assert classify({"insertion_right_edge": True, "reference_DQA2": True}) == "not_classified"
    assert classify({"reference_deletion_point": True, "reference_DQA2": True,
                     "reference_DQB_DQA5": True}) == "reference"
