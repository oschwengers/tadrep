import tadrep.config as cfg
import tadrep.io as tio


def extract():
    plasmids = []
    for file in cfg.plasmids_to_extract:
        circular_genomes = tio.import_sequences(file)
        for id, genome in circular_genomes.items():
            if('circular=true' not in genome['description']):
                circular_genomes.pop(id)
                continue
            genome['file'] = file.name
            
        longest_genome = sorted(circular_genomes.values(), key=lambda k: k['length'])[-1]

        if(longest_genome['length'] > 112_000):
            circular_genomes.pop(longest_genome['id'])
        plasmids.append(circular_genomes)
    
    return plasmids