import logging
import tadrep.config as cfg
from Bio.Seq import Seq

log = logging.getLogger('PLASMIDS')


def detect_plasmids(filtered_contigs, plasmids):
    detected_plasmids = []

    for plasmid_id, hits in filtered_contigs.items():
        plasmid = plasmids[plasmid_id]
        plasmid_coverage, plasmid_identity = calc_coverage_identity(plasmid, hits)

        if (plasmid_coverage >= cfg.min_plasmid_coverage):
            if (plasmid_identity >= cfg.min_plasmid_identity):
                detected_plasmid = {
                    'id': plasmid_id,
                    'hits': hits,
                    'coverage': plasmid_coverage,
                    'identity': plasmid_identity
                }
                detected_plasmids.append(detected_plasmid)
                log.info("plasmid detected: id=%s, contig-hits=%s, total identity=%s, total coverage=%s", plasmid_id, len(hits), plasmid_identity, plasmid_coverage)
            else:
                log.debug("plasmid mismatch: plasmid-id=%s, identity=%s", plasmid_id, plasmid_identity)
        else:
            log.debug("plasmid mismatch: plasmid-id=%s, coverage=%s", plasmid_id, plasmid_coverage)
    return detected_plasmids


def reconstruct_plasmid(plasmid, genome, contigs):
    log.debug('reconstruct plasmid: id=%s, contigs=%s', plasmid['id'], len(plasmid['hits']))

    plasmid['hits'] = [contig for contig in sorted(plasmid['hits'], key=lambda k: k['plasmid_start'])]
    matched_contigs_sorted = [contigs[hit['contig_id']] for hit in plasmid['hits']]
    log.info('sorted contigs: total=%s', len(matched_contigs_sorted))

    sequences = []
    for count, contig in enumerate(plasmid['hits']):
        if(contig['strand'] == '+'):
            sequences.append(matched_contigs_sorted[count]['sequence'])
        else:
            sequences.append(str(Seq(matched_contigs_sorted[count]['sequence']).reverse_complement()))

    gap_sequence = 'N' * cfg.gap_sequence_length
    sequence = gap_sequence.join(sequences)

    plasmid['sequence'] = sequence
    plasmid['description'] = f"{genome}_{plasmid['id']} pseudo plasmid reference={plasmid['id']} contigs={len(plasmid['hits'])}, coverage={plasmid['coverage']:.3f}, identity={plasmid['identity']:.3f}"
    log.info('reconstruct plasmid: id=%s, length=%s, description=%s', plasmid['id'], len(plasmid['sequence']), plasmid['description'])
    return matched_contigs_sorted


def calc_coverage_identity(plasmid, hits):

    contig_length_sum = sum([contig['length'] for contig in hits])
    plasmid_coverage = contig_length_sum / plasmid['length']

    sum_contig_identity = sum([contig['num_identity'] for contig in hits])
    sum_contig_alignment = sum([contig['length'] for contig in hits])
    plasmid_identity = sum_contig_identity / sum_contig_alignment

    return plasmid_coverage, plasmid_identity
