import logging

from Bio.Seq import Seq

import tadrep.config as cfg


log = logging.getLogger('PLASMIDS')


def detect_reference_plasmids(sample, filtered_hits, reference_plasmids):
    detected_plasmids = []
    for reference_plasmid_id, hits in filtered_hits.items():
        reference_plasmid = reference_plasmids[reference_plasmid_id]
        coverage, covered_bp, uncovered_bp = calc_coverage(reference_plasmid, hits)
        identity = calc_identity(hits)
        log.debug("reference plasmid hit: id=%s, identity=%.3f, coverage=%.3f, covered=%i bp, uncovered=%i bp", reference_plasmid_id, identity, coverage, covered_bp, uncovered_bp)
        if (coverage >= cfg.min_plasmid_coverage and identity >= cfg.min_plasmid_identity):
            detected_plasmid = {
                'id': f"{sample}_{reference_plasmid_id}",
                'reference': reference_plasmid_id,
                'genome': sample,
                'hits': hits,
                'coverage': coverage,
                'covered_bp': covered_bp,
                'uncovered_bp': uncovered_bp,
                'identity': identity,
                'length': reference_plasmid['length']
            }
            detected_plasmids.append(detected_plasmid)
            log.info("plasmid detected: id=%s, contig-hits=%i, total identity=%.3f, total coverage=%.3f", reference_plasmid_id, len(hits), identity, coverage)
    return detected_plasmids


def reconstruct_plasmid(plasmid, contigs):
    log.debug('reconstruct plasmid: genome=%s, id=%s, # contigs=%i', plasmid['genome'], plasmid['id'], len(plasmid['hits']))
    plasmid['hits'] = sorted(plasmid['hits'], key=lambda k: k['reference_plasmid_start'])
    sorted_plasmid_contigs = []
    sorted_plasmid_contig_sequences = []
    for hit in plasmid['hits']:
        contig = contigs[hit['contig_id']]
        sorted_plasmid_contigs.append(contig)
        sorted_plasmid_contig_sequences.append(contig['sequence'] if hit['strand'] == '+' else str(Seq(contig['sequence']).reverse_complement()))

    gap_sequence = 'N' * cfg.gap_sequence_length
    plasmid['sequence'] = gap_sequence.join(sorted_plasmid_contig_sequences)
    plasmid['description'] = f"reference={plasmid['id']} contigs={len(plasmid['hits'])} coverage={plasmid['coverage']:.3f} identity={plasmid['identity']:.3f}"
    log.info('plasmid reconstructed: genome=%s, id=%s, length=%s, description=%s', plasmid['genome'], plasmid['id'], len(plasmid['sequence']), plasmid['description'])
    return sorted_plasmid_contigs


def calc_coverage(plasmid, hits):
    cov_array = [0 for i in range(0,plasmid['length'])]
    for hit in hits:
        for i in range(hit['reference_plasmid_start']-1, hit['reference_plasmid_end']):
            cov_array[i] += 1
    uncovered_bp = cov_array.count(0)
    covered_bp = plasmid['length'] - uncovered_bp
    coverage = covered_bp / plasmid['length']
    log.debug("coverage: cov=%.3f, covered-bp=%i, uncovered-bp=%i", coverage, covered_bp, uncovered_bp)
    return coverage, covered_bp, uncovered_bp


def calc_identity(hits):
    sum_num_identity = sum([hit['num_identity'] for hit in hits])
    sum_hit_alignment = sum([hit['length'] for hit in hits])
    identity = sum_num_identity / sum_hit_alignment
    return identity
