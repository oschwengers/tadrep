import logging

import tadrep.config as cfg
import tadrep.utils as tu


log = logging.getLogger('BLAST')


############################################################################
# Setup and run blastn search
############################################################################
def search_contigs(genome_path, blast_output_path):
    db_target = '-subject' if cfg.plasmids_path else '-db'
    db_target_path = cfg.plasmids_path if cfg.plasmids_path else cfg.database_path.joinpath('db')

    cmd_blast = [
        'blastn',
        '-query', str(genome_path),
        db_target, str(db_target_path),
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
            cols = line.strip().split('\t')
            log.debug('raw hit: %s', cols)
            if(cfg.database_path):
                if('|' in cols[4]):
                    cols[4] = cols[4].split('|')[1]
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
                'num_identity': int(cols[8]),
                'evalue': float(cols[10]),
                'bitscore': float(cols[11])
            }
            hits.append(hit)
    log.info('raw blast hits: genome=%s, # hits=%i', genome_path.stem, len(hits))
    return hits


############################################################################
# Parse and filter contig hits
############################################################################
def filter_contig_hits(raw_hits):
    filtered_hits = {}
    for hit in raw_hits:
        if(hit['strand'] == '-'):
            hit['reference_plasmid_start'], hit['reference_plasmid_end'] = hit['reference_plasmid_end'], hit['reference_plasmid_start']

        if(hit['coverage'] >= cfg.min_contig_coverage and hit['perc_identity'] >= cfg.min_contig_identity):
            plasmid_hits = filtered_hits.get(hit['reference_plasmid_id'], [])
            plasmid_hits.append(hit)
            if(len(plasmid_hits) == 1):
                filtered_hits[hit['reference_plasmid_id']] = plasmid_hits
            log.debug(
                'filtered hit: contig-id=%s, reference-plasmid-id=%s, alignment-length=%i, identity=%0.3f, contig-coverage=%0.3f',
                hit['contig_id'], hit['reference_plasmid_id'], hit['length'], hit['perc_identity'], hit['coverage']
            )

    log.info('filtered blast hits: genome=%s, # raw-hits=%i, # filtered-hits=%i', len(raw_hits), len(filtered_hits))
    return filtered_hits
