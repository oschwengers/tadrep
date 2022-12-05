import logging
import sys
import re

import pyrodigal

import tadrep.io as tio
import tadrep.config as cfg
import tadrep.utils as tu


log = logging.getLogger('CHARACTERIZE')


def characterize():
    # load json from output path
    db_path = cfg.output_path.joinpath('db.json')
    existing_data = tio.load_data(db_path)
    plasmids = existing_data.get('plasmids', {})

    # check if data is available
    if(not existing_data or not plasmids):
        log.error('Failed to load data from %s', db_path)
        sys.exit(f'ERROR: Failed to load data from {db_path}!')

    # write multifasta
    fasta_path = cfg.output_path.joinpath('db.fasta')
    tio.export_sequences(plasmids.values(), fasta_path)

    # search inc_types for all plamids
    download_inc_reference()
    inc_types_per_plasmid = search_inc_types(fasta_path)

    # write multifasta
    fasta_path = cfg.output_path.joinpath('db.fasta')
    tio.export_sequences(plasmids.values(), fasta_path)

    search_inc_types(fasta_path)

    for plasmid in plasmids.values():
        # set length
        plasmid['length'] = len(plasmid['sequence'])

        # set gc_content
        plasmid['gc_content'] = gc_content(plasmid['sequence'])

        # set individual INC_types
        plasmid['inc_types'] = inc_types_per_plasmid.get(plasmid['id'], [])

        # gene prediction (pyrodigal)
        plasmid['cds'] = gene_prediction(plasmid['sequence'])

        cfg.verboseprint(f"Plasmid: {plasmid['id']:20} Length: {plasmid['length']:7} GC: {plasmid['gc_content']:.2} CDS: {len(plasmid['cds']):5} INC_Types: {len(plasmid['inc_types']):3}")
        log.info('Plasmid: %s, len: %d, gc: %f, cds: %d, inc_types: %d', plasmid['id'], plasmid['length'], plasmid['gc_content'], len(plasmid['cds']), len(plasmid['inc_types']))

    download_inc_types()
    # update json
    tio.export_json(existing_data, db_path)


def gc_content(sequence):
    # count C and G occurenxes
    c_content = sequence.count('C')
    g_content = sequence.count('G')
    gc_combined = c_content + g_content

    # calculate gc percentage
    percentage = gc_combined / float(len(sequence))

    return percentage


def download_inc_reference():

    inc_types_path = cfg.output_path.joinpath('inc-types.fasta')

    # check if inc-types.fasta is already available
    if(inc_types_path.is_file()):
        cfg.verboseprint(f'Found {inc_types_path.name} in {cfg.output_path}')
        log.debug('Found inc-types.fasta in path %s', cfg.output_path)
        return

    plasmidfinder_git = 'https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git'

    cfg.verboseprint(f'Missing {inc_types_path.name} in output directory, attempting download... ')
    log.debug('Missing %s, cloning from %s', inc_types_path.name, plasmidfinder_git)
    
    git_cmd = ['git', 'clone', str(plasmidfinder_git)]
    tu.run_cmd(git_cmd, cfg.tmp_path)

    with open(inc_types_path, 'w') as inc_types:

        # read each .fsa file
        plasmidfinder_path = cfg.tmp_path.joinpath('plasmidfinder_db')
        for fasta_file in plasmidfinder_path.glob('*.fsa'):
            log.debug('Reading %s', fasta_file)
            with open(fasta_file, 'r') as inc_raw:

                for line in inc_raw.readlines():
                    # if line is header
                    if(line.startswith('>')):
                        # check if pattern
                        if(re.search('^>([a-zA-Z0-9()]+)_([0-9]+)_(.*)_(.*)$', line)):
                            inc_types.write(line.split('_')[0])
                        else:
                            inc_types.write(line)
                    # write seq in uppercase
                    else:
                        inc_types.write(line.upper())
                    inc_types.write('\n')


def search_inc_types(db_path):
    """Search for incompatibility motifs."""

    tmp_output_path = cfg.output_path.joinpath('db.inc.blast.out')

    inc_types_cmd = [
        'blastn',
        '-query', str(cfg.output_path.joinpath('inc-types.fasta')),
        '-subject', str(db_path),
        '-num_threads', str(cfg.threads),
        '-perc_identity', '90',
        '-culling_limit', '1',
        '-outfmt', '6 qseqid sseqid sstart send sstrand pident qcovs bitscore',
        '-out', str(tmp_output_path)
    ]

    tu.run_cmd(inc_types_cmd, cfg.output_path)

    hits_per_plasmid = {}
    with tmp_output_path.open() as fh:
        for line in fh:
            cols = line.rstrip().split('\t')
            hit = {
                'type': cols[0],
                'start': int(cols[2]),
                'end': int(cols[3]),
                'strand': '+' if cols[4] == 'plus' else '-',
                'identity': float(cols[5]) / 100,
                'coverage': float(cols[6]) / 100,
                'bitscore': int(cols[7])
            }
            plasmid_id = cols[1]
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

    return hits_per_plasmid


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


def download_inc_types():

    inc_types_path = cfg.output_path.joinpath('inc-types.fasta')

    # check if inc-types.fasta is already available
    if(inc_types_path.is_file()):
        return

    return
