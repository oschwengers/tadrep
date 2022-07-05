import logging

import tadrep.config as cfg
import tadrep.io as tio

log = logging.getLogger('EXTRACT')


def extract():
    plasmids = []
    for file in cfg.plasmids_to_extract:
        genomes = tio.import_sequences(file)
        not_circular_ids = []
        for id, genome in genomes.items():
            if('circular=true' not in genome['description']):
                not_circular_ids.append(id)
                continue
            genome['file'] = file.name
        
        log.info('loaded genomes: total=%i, not circular=%i file=%s', len(genomes), len(not_circular_ids), file.name)
        circular_genomes = {k: genomes[k] for k in genomes.keys() if k not in not_circular_ids}

        if(not circular_genomes):
            continue

        longest_genome = sorted(circular_genomes.values(), key=lambda k: k['length'])[-1]
        log.info('longest genome: id=%s, length=%i', longest_genome['id'], longest_genome['length'])

        if(longest_genome['length'] > 112_000):
            circular_genomes.pop(longest_genome['id'])
            log.debug('dropped genome=%s, length=%i', longest_genome['id'], longest_genome['length'])
        plasmids.append(circular_genomes)
    
    return plasmids