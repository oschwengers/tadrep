from pathlib import Path
from unittest.mock import patch

import tadrep.blast as tb


@patch("tadrep.blast.cfg.min_contig_coverage", 40 / 100)
@patch("tadrep.blast.cfg.min_contig_identity", 90 / 100)
def test_hit_filtering():
    raw_hits = []
    with Path("test/data/blastn.tsv").open('r') as fh:
        for line in fh:
            cols = line.strip().split('\t')
            hit = {
                'contig_id': cols[0],
                'contig_start': int(cols[1]),
                'contig_end': int(cols[2]),
                'contig_length': int(cols[3]),
                'reference_plasmid_id': cols[4],
                'reference_plasmid_start': int(cols[5]),
                'reference_plasmid_end': int(cols[6]),
                'length': int(cols[7]),
                'strand': '+' if cols[9] == 'plus' else '-',
                'coverage': int(cols[7]) / int(cols[3]),
                'perc_identity': int(cols[8]) / int(cols[7]),
                'num_identity': int(cols[8])
            }
            raw_hits.append(hit)
    assert len(raw_hits) == 36

    expected_hits = {
        "p1": [
            15, 23, 27, 31, 35, 36, 40, 48, 51, 53, 57, 58, 67, 68, 78
        ],
        "p2": [
            24, 32, 41, 49, 60
        ],
        "p3": [
            34, 34
        ]
    }

    # assert with hits
    filtered_hits = tb.filter_contig_hits('test', raw_hits)
    assert filtered_hits.keys() == expected_hits.keys()
    assert len(filtered_hits) == len(expected_hits)

    for plasmid_id, hits in expected_hits.items():
        assert len(filtered_hits[plasmid_id]) == len(hits)

    # assert without hits
    filtered_hits = tb.filter_contig_hits('test', {})
    assert filtered_hits == {}
