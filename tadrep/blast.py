import logging

import tadrep.config as cfg
import tadrep.utils as tu


log = logging.getLogger('BLAST')


############################################################################
# Setup and run blastn search
############################################################################
def search_contigs(genome_path, blast_output_path):

    cmd_blast = [
        'blastn',
        '-query', str(genome_path),
        '-db', str(cfg.output_path.joinpath('reference_db/db')),
        '-culling_limit', '1',
        '-evalue', '1E-5',
        '-num_threads', str(cfg.blast_threads),
        '-outfmt', '6 qseqid qstart qend qlen sseqid sstart send length nident sstrand evalue bitscore',
        '-out', str(blast_output_path)
    ]
    log.debug('cmd=%s', cmd_blast)

    tu.run_cmd(cmd_blast, cfg.tmp_path)

    hits = []
    with blast_output_path.open('r') as fh:
        for line in fh:
            (qseqid, qstart, qend, qlen, sseqid, sstart, send, length, nident, sstrand, evalue, bitscore) = line.strip().split('\t')
            hit = {
                'contig_id': qseqid,
                'contig_start': int(qstart),
                'contig_end': int(qend),
                'contig_length': int(qlen),
                'reference_plasmid_id': sseqid,
                'reference_plasmid_start': int(sstart),
                'reference_plasmid_end': int(send),
                'length': int(length),
                'strand': '+' if sstrand == 'plus' else '-',
                'coverage': int(length) / int(qlen),
                'perc_identity': int(nident) / int(length),
                'num_identity': int(nident),
                'evalue': float(evalue),
                'bitscore': float(bitscore)
            }
            if(hit['strand'] == '-'):
                hit['reference_plasmid_start'], hit['reference_plasmid_end'] = hit['reference_plasmid_end'], hit['reference_plasmid_start']
            hits.append(hit)
    log.info('raw blast hits: genome=%s, # hits=%i', genome_path.stem, len(hits))
    return hits


############################################################################
# Parse and filter contig hits
############################################################################
def filter_contig_hits(genome, raw_hits, reference_plasmids):
    filtered_hits = 0
    filtered_hits_per_ref_plasmid = {}
    edge_hits_per_ref_plasmid = {}
    for hit in raw_hits:
        reference_plasmid_id = hit['reference_plasmid_id']
        if(hit['perc_identity'] >= cfg.min_contig_identity):
            reference_plasmid = reference_plasmids[reference_plasmid_id]
            if(hit['reference_plasmid_start'] == 1 or hit['reference_plasmid_end'] == reference_plasmid['length']):  # hit at plasmid edge either 5' or 3', store for combined coverage check
                edge_hits_per_contig = edge_hits_per_ref_plasmid.get(reference_plasmid_id, [])
                edge_hits_per_contig.append(hit)
                if(reference_plasmid_id not in edge_hits_per_ref_plasmid):
                    edge_hits_per_ref_plasmid[reference_plasmid_id] = edge_hits_per_contig
            elif(hit['coverage'] >= cfg.min_contig_coverage):  # hit within plasmid with sufficient coverage
                plasmid_hits = filtered_hits_per_ref_plasmid.get(reference_plasmid_id, [])
                plasmid_hits.append(hit)
                filtered_hits += 1
                if(reference_plasmid_id not in filtered_hits_per_ref_plasmid):
                    filtered_hits_per_ref_plasmid[reference_plasmid_id] = plasmid_hits
                log.debug(
                    'filtered hit: contig-id=%s, reference-plasmid-id=%s, alignment-length=%i, identity=%0.3f, contig-coverage=%0.3f',
                    hit['contig_id'], reference_plasmid_id, hit['length'], hit['perc_identity'], hit['coverage']
                )
    
    for reference_plasmid_id, edge_hits in edge_hits_per_ref_plasmid.items():
        if(len(edge_hits) == 1):
            edge_hit = edge_hits[0]
            if(edge_hit['coverage'] >= cfg.min_contig_coverage):
                plasmid_hits = filtered_hits_per_ref_plasmid.get(reference_plasmid_id, [])
                plasmid_hits.append(edge_hit)
                filtered_hits += 1
                if(reference_plasmid_id not in filtered_hits_per_ref_plasmid):
                    filtered_hits_per_ref_plasmid[reference_plasmid_id] = plasmid_hits
                log.debug(
                    'filtered single edge hit: contig-id=%s, reference-plasmid-id=%s, length=%i, identity=%0.3f, contig-coverage=%0.3f, plasmid-start=%i, plasmid-end=%i',
                    edge_hit['contig_id'], reference_plasmid_id, edge_hit['length'], edge_hit['perc_identity'], edge_hit['coverage'], edge_hit['reference_plasmid_start'], edge_hit['reference_plasmid_end']
                )
        elif(len(edge_hits) == 2):
            (edge_hit_a, edge_hit_b) = edge_hits
            if(edge_hit_a['contig_id'] == edge_hit_b['contig_id']):  # check hits belong to the same contig
                alignment_sum = edge_hit_a['length'] + edge_hit_b['length']
                reference_plasmid = reference_plasmids[reference_plasmid_id]
                contig_cov = alignment_sum / edge_hit_a['contig_length']
                contig_ident = (edge_hit_a['num_identity'] + edge_hit_b['num_identity']) / ((edge_hit_a['length'] + edge_hit_b['length']))
                if(contig_cov >= cfg.min_contig_coverage):
                    plasmid_hits = filtered_hits_per_ref_plasmid.get(reference_plasmid_id, [])
                    plasmid_hits.append(edge_hit_a)
                    plasmid_hits.append(edge_hit_b)
                    filtered_hits += 2
                    if(reference_plasmid_id not in filtered_hits_per_ref_plasmid):
                        filtered_hits_per_ref_plasmid[reference_plasmid_id] = plasmid_hits
                    log.debug(
                        'filtered combined edge hits: contig-id=%s, reference-plasmid-id=%s, combined-length=%i, identity=%0.3f, combined-coverage=%0.3f',
                        edge_hit_a['contig_id'], edge_hit_a['reference_plasmid_id'], alignment_sum, contig_ident, contig_cov
                    )

    log.info('filtered blast hits: genome=%s, # raw-hits=%i, # filtered-hits=%i', genome, len(raw_hits), filtered_hits)
    return filtered_hits_per_ref_plasmid
