import logging

from Bio.Seq import Seq

import tadrep.config as cfg


log = logging.getLogger('PLASMIDS')


def detect_plasmids(filtered_contigs, plasmids):
    detected_plasmids = []

    for plasmid_id, hits in filtered_contigs.items():
        plasmid = plasmids[plasmid_id]
        plasmid_coverage, covered_bp, uncovered_bp = calc_coverage(plasmid, hits)
        plasmid_identity = calc_identity(hits)

        if (plasmid_coverage >= cfg.min_plasmid_coverage):
            if (plasmid_identity >= cfg.min_plasmid_identity):
                detected_plasmid = {
                    'id': plasmid_id,
                    'reference': plasmid_id,
                    'hits': hits,
                    'coverage': plasmid_coverage,
                    'covered_bp': covered_bp,
                    'uncovered_bp': uncovered_bp,
                    'identity': plasmid_identity,
                    'length': plasmid['length']
                }
                detected_plasmids.append(detected_plasmid)
                log.info("plasmid detected: id=%s, contig-hits=%s, total identity=%s, total coverage=%s", plasmid_id, len(hits), plasmid_identity, plasmid_coverage)
            else:
                log.debug("plasmid mismatch: plasmid-id=%s, identity=%s", plasmid_id, plasmid_identity)
        else:
            log.debug("plasmid mismatch: plasmid-id=%s, coverage=%s", plasmid_id, plasmid_coverage)
    return detected_plasmids


def reconstruct_plasmid(plasmid, genome, contigs):
    log.debug('reconstruct plasmid: genome=%s, id=%s, contigs=%s', genome, plasmid['id'], len(plasmid['hits']))

    plasmid['hits'] = [contig for contig in sorted(plasmid['hits'], key=lambda k: k['plasmid_start'])]
    plasmid_contigs_sorted = [contigs[hit['contig_id']] for hit in plasmid['hits']]
    log.debug('sorted contigs: total=%s', len(plasmid_contigs_sorted))

    sequences = []
    for count, contig in enumerate(plasmid['hits']):
        if(contig['strand'] == '+'):
            contig_seq = plasmid_contigs_sorted[count]['sequence']
        else:
            contig_seq = str(Seq(plasmid_contigs_sorted[count]['sequence']).reverse_complement())
        sequences.append(contig_seq)

    gap_sequence = 'N' * cfg.gap_sequence_length
    plasmid['sequence'] = gap_sequence.join(sequences)

    plasmid['id'] = f"{genome}_{plasmid['id']}"
    plasmid['description'] = f"reference={plasmid['id']} contigs={len(plasmid['hits'])} coverage={plasmid['coverage']:.3f} identity={plasmid['identity']:.3f}"
    log.info('reconstruct plasmid: genome=%s, id=%s, length=%s, description=%s', genome, plasmid['id'], len(plasmid['sequence']), plasmid['description'])
    return plasmid_contigs_sorted


def calc_coverage(plasmid, hits):
    cov_array = [0 for i in range(0,plasmid['length'])]
    for hit in hits:
        for i in range(hit['plasmid_start']-1, hit['plasmid_end']):
            cov_array[i] += 1
    uncovered_bp = cov_array.count(0)
    covered_bp = plasmid['length'] - uncovered_bp
    coverage = covered_bp / plasmid['length']
    log.debug("coverage: cov=%.3f, covered-bp=%i, uncovered-bp=%i", coverage, covered_bp, uncovered_bp)
    return coverage, covered_bp, uncovered_bp


def calc_identity(hits):
    sum_hit_identity = sum([hit['num_identity'] for hit in hits])
    sum_hit_alignment = sum([hit['length'] for hit in hits])
    plasmid_identity = sum_hit_identity / sum_hit_alignment
    return plasmid_identity
