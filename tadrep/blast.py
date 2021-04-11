import subprocess as sp
import logging
from pathlib import Path

import tadrep.config as cfg

log = logging.getLogger('BLAST')


############################################################################
# Setup blastn search
# Run Blastn search
############################################################################
def search_contigs(genome_path, plasmid_path):
    blast_output_path = cfg.output_path.joinpath('blastn.tsv')

    cmd_blast = [
        'blastn',
        '-query', str(genome_path),
        '-subject', str(plasmid_path),
        '-culling_limit', '1',
        '-evalue', '1E-5',
        '-num_threads', str(cfg.threads),
        '-outfmt', '6 qseqid qstart qend qlen sseqid sstart send length nident sstrand evalue bitscore',
        '-out', str(blast_output_path)
    ]
    log.debug("cmd=%s", cmd_blast)

    process = sp.run(
        cmd_blast,
        cwd=str(cfg.tmp_path),
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )

    if(process.returncode != 0):
        log.debug('stdout=%s, stderr=%s', process.stdout, process.stderr)
        log.warning('blastn failed! Blastn-error-code: %s', process.returncode)
        raise Exception(f'blastn error! Error code: {process.returncode}')

    hits = []
    with blast_output_path.open('r') as fh:
        for line in fh:
            cols = line.strip().split('\t')
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

    log.debug("Total contigs matched: %s", len(hits))

    return hits

############################################################################
# Process contig-hits into columns
# Check strand direction
# return hits
############################################################################
def parse_contig_hits():
    pass
