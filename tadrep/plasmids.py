import logging

from Bio.Seq import Seq

import tadrep.config as cfg


log = logging.getLogger('PLASMIDS')


def detect_plasmids(filtered_contigs, plasmids):
    detected_plasmids = []

    for plasmid_id, hits in filtered_contigs.items():
        plasmid = plasmids[plasmid_id]
        plasmid_coverage = calc_coverage(plasmid, hits)
        plasmid_identity = calc_identity(hits)

        if (plasmid_coverage >= cfg.min_plasmid_coverage):
            if (plasmid_identity >= cfg.min_plasmid_identity):
                detected_plasmid = {
                    'id': plasmid_id,
                    'reference': plasmid_id,
                    'hits': hits,
                    'coverage': plasmid_coverage,
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
    hits_length_sum = sum([hit['length'] for hit in hits])
    plasmid_coverage = hits_length_sum / plasmid['length']
    return plasmid_coverage


def calc_identity(hits):
    sum_hit_identity = sum([hit['num_identity'] for hit in hits])
    sum_hit_alignment = sum([hit['length'] for hit in hits])
    plasmid_identity = sum_hit_identity / sum_hit_alignment
    return plasmid_identity
