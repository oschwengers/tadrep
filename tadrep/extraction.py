import tadrep.config as cfg
import tadrep.io as tio


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
        
        circular_genomes = {k: genomes[k] for k in genomes.keys() if k not in not_circular_ids}

        if(not circular_genomes):
            continue

        longest_genome = sorted(circular_genomes.values(), key=lambda k: k['length'])[-1]

        if(longest_genome['length'] > 112_000):
            circular_genomes.pop(longest_genome['id'])
        plasmids.append(circular_genomes)
    
    return plasmids