from unittest.mock import patch

import tadrep.plasmids as tp


@patch('tadrep.plasmids.cfg.min_plasmid_coverage', 80 / 100)
@patch('tadrep.plasmids.cfg.min_plasmid_identity', 90 / 100)
def test_plasmid_detection():
    plasmids = {
        'test_p1': {
            'id': 'test_p1',
            'length': 5000
        },
        'test_p2': {
            'id': 'test_p2',
            'length': 5000
        },
        'test_p3': {
            'id': 'test_p3',
            'length': 5400
        }
    }
    filtered_contigs = {
        # Identity and Coverage above threshold
        'test_p1': [
            {
                'contig_id': 'test_c1',
                'length': 1300,
                'num_identity': 1200
            },
            {
                'contig_id': 'test_c2',
                'length': 2100,
                'num_identity': 1900
            },
            {
                'contig_id': 'test_c3',
                'length': 900,
                'num_identity': 899
            },
        ],
        # Identity below threshold
        'test_p2': [
            {
                'contig_id': 'test_c1',
                'length': 1300,
                'num_identity': 1100
            },
            {
                'contig_id': 'test_c2',
                'length': 2100,
                'num_identity': 1800
            },
            {
                'contig_id': 'test_c3',
                'length': 900,
                'num_identity': 700
            },
        ],
        # Coverage below threshold
        'test_p3': [
            {
                'contig_id': 'test_c1',
                'length': 1300,
                'num_identity': 1300
            },
            {
                'contig_id': 'test_c2',
                'length': 2100,
                'num_identity': 2100
            },
            {
                'contig_id': 'test_c3',
                'length': 900,
                'num_identity': 900
            },
        ]
    }
    expected_plasmids = [
        {
            'id': 'test_p1',
            'reference': 'test_p1',
            'hits': filtered_contigs['test_p1'],
            'coverage': 0.86,
            'identity': 0.93,
            'length': 5000
        }
    ]
    detected_plasmids = tp.detect_plasmids(filtered_contigs, plasmids)
    assert detected_plasmids == expected_plasmids

    detected_plasmids = tp.detect_plasmids({}, plasmids)
    assert detected_plasmids == []
