import tadrep.plasmids as tp
from unittest.mock import patch


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
            'length': 2000
        }
    }
    filtered_contigs = {
        'test_p1': [
            {
                'contig_id': 'test_c1',
                'length': 100,
                'num_identity': 90
            },
            {
                'contig_id': 'test_c2',
                'length': 100,
                'num_identity': 90
            },
            {
                'contig_id': 'test_c3',
                'length': 100,
                'num_identity': 90
            },
        ],
        'test_p2': [
            {
                'contig_id': 'test_c1',
                'length': 100,
                'num_identity': 90
            },
            {
                'contig_id': 'test_c1',
                'length': 100,
                'num_identity': 90
            },
            {
                'contig_id': 'test_c1',
                'length': 100,
                'num_identity': 90
            },
        ],
        'test_p3': []
    }
