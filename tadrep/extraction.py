import logging

import tadrep.config as cfg
import tadrep.io as tio

log = logging.getLogger('EXTRACT')


def extract():
    plasmids = {}
    for file in cfg.plasmids_to_extract:
        genomes = tio.import_sequences(file, sequence=True)
        
        circular_genomes = get_circular(genomes)

        if(not circular_genomes):
            continue

        if(cfg.drop_longest):
            circular_genomes = filter_longest(circular_genomes)

        circular_genomes = update_keys(circular_genomes, file.name)
        plasmids.update(circular_genomes)
    
    plasmid_output_path = cfg.output_path.joinpath('plasmids.json')
    tio.export_json(plasmids, plasmid_output_path)


def get_circular(genomes):
    not_circular_ids = []
    for id, genome in genomes.items():
        if('circular=true' not in genome['description']):
            not_circular_ids.append(id)
            continue
    
    log.info('loaded genomes: total=%i, not circular=%i', len(genomes), len(not_circular_ids))
    circular_genomes = {k: genomes[k] for k in genomes.keys() if k not in not_circular_ids}

    return circular_genomes


def filter_longest(circular_genomes):
        longest_genome = sorted(circular_genomes.values(), key=lambda k: k['length'])[-1]
        log.info('longest genome: id=%s, length=%i', longest_genome['id'], longest_genome['length'])

        if(longest_genome['length'] > 112_000):
            circular_genomes.pop(longest_genome['id'])
            log.debug('dropped genome=%s, length=%i', longest_genome['id'], longest_genome['length'])
        return circular_genomes


def update_keys(circular_dict, filename):
    for key in circular_dict.keys():
        circular_dict[key]['filename'] = filename

    updated_dict = dict((f'{filename}_p{plasmid_id}', value) for plasmid_id, value in circular_dict.items())

    return updated_dict
