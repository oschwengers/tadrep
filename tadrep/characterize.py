import logging
import sys

import pyrodigal

import tadrep.io as tio
import tadrep.config as cfg
import tadrep.utils as tu


log = logging.getLogger('CHARACTERIZE')


def characterize():
    # load json from output path
    db_path = cfg.output_path.joinpath('db.json')
    db_data = tio.load_data(db_path)

    # check if data is available
    if(not db_data):
        log.error('Failed to load data from %s', db_path)
        sys.exit(f'ERROR: Failed to load data from {db_path}!')
    if(not db_data.get('plasmids', {})):
        log.error('Failed to load Plasmids from %s', db_path)
        sys.exit(f'ERROR: Failed to load Plasmids from {db_path}! Maybe file is empty?')

    # write multifasta
    fasta_path = cfg.tmp_path.joinpath('db.fasta')
    tio.export_sequences(db_data['plasmids'].values(), fasta_path)

    # search inc_types for all plasmids
    inc_types_per_plasmid = search_inc_types(fasta_path)

    for id, plasmid in db_data['plasmids'].items():
        plasmid['length'] = len(plasmid['sequence'])
        plasmid['gc_content'] = calc_gc_content(plasmid['sequence'])
        plasmid['inc_types'] = inc_types_per_plasmid.get(plasmid['id'], [])
        plasmid['cds'] = gene_prediction(plasmid['sequence'])  # gene prediction

        inc_types = ','.join([inc['type'] for inc in plasmid['inc_types']])
        cfg.verbose_print(f"{plasmid['id']}:\tLength: {plasmid['length']:8}\tGC: {plasmid['gc_content']:3.2}\tCDS: {len(plasmid['cds']):5}\tInc Types: {inc_types}")
        log.info('Plasmid: %s, len: %d, gc: %f, cds: %d, inc_types: %d', plasmid['id'], plasmid['length'], plasmid['gc_content'], len(plasmid['cds']), len(plasmid['inc_types']))

    # update json
    print('Writing JSON...')
    tio.export_json(db_data, db_path)


def calc_gc_content(sequence):
    return (sequence.count('C') + sequence.count('G')) / float(len(sequence))


def search_inc_types(db_path):
    """Search for incompatibility motifs."""
    inc_types = cfg.output_path.joinpath('inc-types.fasta')
    if(not inc_types.is_file()):
        log.debug("Inc_types reference not found!")
        sys.exit("ERROR: Inc_types reference not found! Please import with '--inc-types PATH_TO_FASTA' or download it with subcommand setup!")

    tmp_output_path = cfg.tmp_path.joinpath('db.inc.blast.out')
    inc_types_cmd = [
        'blastn',
        '-query', str(inc_types),
        '-subject', str(db_path),
        '-num_threads', str(cfg.threads),
        '-perc_identity', '90',
        '-culling_limit', '1',
        '-outfmt', '6 qseqid sseqid sstart send sstrand pident qcovs bitscore',
        '-out', str(tmp_output_path)
    ]
    tu.run_cmd(inc_types_cmd, cfg.output_path)

    hits_per_plasmid = {}
    with tmp_output_path.open('r') as fh:
        for line in fh:
            cols = line.rstrip().split('\t')
            plasmid_id = cols[1]
            hit = {
                'type': cols[0],
                'start': int(cols[2]),
                'end': int(cols[3]),
                'strand': '+' if cols[4] == 'plus' else '-',
                'identity': float(cols[5]) / 100,
                'coverage': float(cols[6]) / 100,
                'bitscore': int(float(cols[7]))
            }
            if(hit['coverage'] >= 0.6):
                hits_per_pos = hits_per_plasmid.get(plasmid_id, {})
                hit_pos = hit['end'] if hit['strand'] == '+' else hit['start']
                if(hit_pos in hits_per_pos):
                    former_hit = hits_per_pos[hit_pos]
                    if(hit['bitscore'] > former_hit['bitscore']):
                        hits_per_pos[hit_pos] = hit
                        log.info(
                            'inc-type: hit! contig=%s, type=%s, start=%d, end=%d, strand=%s',
                            plasmid_id, hit['type'], hit['start'], hit['end'], hit['strand']
                        )
                else:
                    hits_per_pos[hit_pos] = hit
                    log.info(
                        'inc-type: hit! contig=%s, type=%s, start=%d, end=%d, strand=%s',
                        plasmid_id, hit['type'], hit['start'], hit['end'], hit['strand']
                    )
                hits_per_plasmid[plasmid_id] = hits_per_pos
    
    filtered_hits_per_plasmid = {}
    for plasmid_id, hits in hits_per_plasmid.items():  # remove potential smaller partial hits
        hits_per_inc = {}
        for hit in hits.values():
            if hit['type'] in hits_per_inc:
                former_hit = hits_per_inc[hit['type']]
                if(hit['bitscore'] > former_hit['bitscore']):
                    hits_per_inc[hit['type']] = hit
            else:
                hits_per_inc[hit['type']] = hit
        filtered_hits_per_plasmid[plasmid_id] = list(hits_per_inc.values())
    
    return filtered_hits_per_plasmid


def gene_prediction(plasmid_sequence):
    plasmid_cds = []
    orffinder = pyrodigal.OrfFinder(meta=True, closed=True)
    for prediction in orffinder.find_genes(plasmid_sequence.encode()):
        hit = {
            'start': prediction.begin,
            'stop': prediction.end,
            'strand': '+' if prediction.strand == 1 else '-',
            'aa': prediction.translate()
        }
        plasmid_cds.append(hit)

    return plasmid_cds
