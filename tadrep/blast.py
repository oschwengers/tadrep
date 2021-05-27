import logging

import tadrep.config as cfg
import tadrep.utils as tu


log = logging.getLogger('BLAST')


############################################################################
# Setup blastn search
# Run Blastn search
############################################################################
def search_contigs(genome_path):
    blast_output_path = cfg.tmp_path.joinpath('blastn.tsv')

    db_target = '-subject' if cfg.plasmids_path else '-db'
    db_target_path = cfg.plasmids_path if cfg.plasmids_path else cfg.database_path.joinpath("db")

    cmd_blast = [
        'blastn',
        '-query', str(genome_path),
        db_target, str(db_target_path),
        '-culling_limit', '1',
        '-evalue', '1E-5',
        '-num_threads', str(cfg.threads),
        '-outfmt', '6 qseqid qstart qend qlen sseqid sstart send length nident sstrand evalue bitscore',
        '-out', str(blast_output_path)
    ]
    log.debug("cmd=%s", cmd_blast)

    tu.run_cmd(cmd_blast, cfg.tmp_path)

    hits = []
    log.debug('reading blast output: file=%s', blast_output_path)
    with blast_output_path.open('r') as fh:
        for line in fh:
            cols = line.strip().split('\t')
            if(cfg.database_path):
                if('|' in cols[4]):
                    cols[4] = cols[4].split('|')[1]
            hit = {
                'contig_id': cols[0],
                'contig_start': int(cols[1]),
                'contig_end': int(cols[2]),
                'contig_length': int(cols[3]),
                'plasmid_id': cols[4],
                'plasmid_start': int(cols[5]),
                'plasmid_end': int(cols[6]),
                'length': int(cols[7]),
                'strand': '+' if cols[9] == 'plus' else '-',
                'coverage': int(cols[7]) / int(cols[3]),
                'perc_identity': int(cols[8]) / int(cols[7]),
                'num_identity': int(cols[8]),
                'evalue': float(cols[10]),
                'bitscore': float(cols[11])
            }
            hits.append(hit)

    log.info('search blast-hits: genome=%s, plasmids=%s, found=%s', genome_path.stem, db_target_path.stem, len(hits))
    return hits


############################################################################
# Process contig-hits into columns
# Check strand direction
# return hits
############################################################################
def filter_contig_hits(raw_hits):
    filtered_hits = {}
    for hit in raw_hits:
        if(hit['strand'] == '-'):
            hit['plasmid_start'], hit['plasmid_end'] = hit['plasmid_end'], hit['plasmid_start']

        if(hit['coverage'] >= cfg.min_contig_coverage and hit['perc_identity'] >= cfg.min_contig_identity):
            plasmid_hits = filtered_hits.get(hit['plasmid_id'], [])
            plasmid_hits.append(hit)
            if(len(plasmid_hits) == 1):
                filtered_hits[hit['plasmid_id']] = plasmid_hits
            log.debug('contig hit: contig-id=%s, plasmid-id=%s, alignment length=%s, identity=%s', hit['contig_id'], hit['plasmid_id'], hit['length'], hit['perc_identity'])

    log.info('filtered blast-hits: raw-hits=%s, filtered-hits=%s', len(raw_hits), len(filtered_hits))
    return filtered_hits
