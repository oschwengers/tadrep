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
                'reference_plasmid_start': 1,
                'reference_plasmid_end': 1300,
                'num_identity': 1200
            },
            {
                'contig_id': 'test_c2',
                'length': 2100,
                'reference_plasmid_start': 1301,
                'reference_plasmid_end': 3400,
                'num_identity': 1900
            },
            {
                'contig_id': 'test_c3',
                'length': 900,
                'reference_plasmid_start': 3401,
                'reference_plasmid_end': 4300,
                'num_identity': 899
            },
        ],
        # Identity below threshold
        'test_p2': [
            {
                'contig_id': 'test_c1',
                'length': 1300,
                'reference_plasmid_start': 1,
                'reference_plasmid_end': 1300,
                'num_identity': 1000
            },
            {
                'contig_id': 'test_c2',
                'length': 2100,
                'reference_plasmid_start': 1301,
                'reference_plasmid_end': 3400,
                'num_identity': 1700
            },
            {
                'contig_id': 'test_c3',
                'length': 900,
                'reference_plasmid_start': 3401,
                'reference_plasmid_end': 4300,
                'num_identity': 500
            },
        ],
        # Coverage below threshold
        'test_p3': [
            {
                'contig_id': 'test_c1',
                'length': 900,
                'reference_plasmid_start': 1,
                'reference_plasmid_end': 900,
                'num_identity': 850
            },
            {
                'contig_id': 'test_c2',
                'length': 1000,
                'reference_plasmid_start': 1301,
                'reference_plasmid_end': 2300,
                'num_identity': 900
            },
            {
                'contig_id': 'test_c3',
                'length': 500,
                'reference_plasmid_start': 3401,
                'reference_plasmid_end': 3900,
                'num_identity': 400
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
    detected_plasmids = tp.detect_reference_plasmids('test-sample', filtered_contigs, plasmids)
    assert len(detected_plasmids) == 1
    detected_plasmid = detected_plasmids[0]
    assert detected_plasmid['coverage'] == expected_plasmids[0]['coverage']
    assert detected_plasmid['identity'] == expected_plasmids[0]['identity']

    detected_plasmids = tp.detect_reference_plasmids('test-sample', {}, plasmids)
    assert detected_plasmids == []
