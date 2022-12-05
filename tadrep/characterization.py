import logging
import sys

import tadrep.io as tio
import tadrep.config as cfg
import tadrep.utils as tu


log = logging.getLogger('CHARACTERIZE')


def characterize():
    # load json from output path
    json_path = cfg.output_path.joinpath('db.json')
    existing_data = tio.load_data(json_path)
    plasmids = existing_data.get('plasmids', {})

    # check if data is available
    if(not existing_data or not plasmids):
        log.error('Failed to load data from %s', json_path)
        sys.exit(f'ERROR: Failed to load data from {json_path}!')

    # write multifasta
    fasta_path = cfg.output_path.joinpath('db.fasta')
    tio.export_sequences(plasmids.values(), fasta_path)

    search_inc_types(fasta_path)

    for plasmid in plasmids.values():
        # set length
        plasmid['length'] = len(plasmid['sequence'])

        # set gc_content
        plasmid['gc_content'] = gc_content(plasmid['sequence'])

        # calculate INC_types (Platon)
            # write multifasta
            # find inc-types.fasta
            # threads = cfg.threads

        # gene prediction (pyrodigal)

    download_inc_types()
    # update json
    tio.export_json(existing_data, json_path)


def gc_content(sequence):
    # count C and G occurenxes
    c_content = sequence.count('C')
    g_content = sequence.count('G')
    gc_combined = c_content + g_content

    # calculate gc percentage
    percentage = gc_combined / float(len(sequence))

    return percentage


def download_inc_types():

    inc_types_path = cfg.output_path.joinpath('inc-types.fasta')

    # check if inc-types.fasta is already available
    if(inc_types_path.is_file()):
        return

    return


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

    return hits_per_plasmid


def gene_prediction():
    pass
