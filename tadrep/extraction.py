import tadrep.config as cfg
import tadrep.io as tio

def extract():
    plasmids = []
    for file in cfg.plasmids_to_extract:
        circular_genomes = tio.import_sequences(file)
        longest_genome = sorted(circular_genomes.values(), key=lambda k: k['length'])[-1]

        if(longest_genome['length'] > 112_000):
            circular_genomes.drop(longest_genome['id'])
        plasmids.append(circular_genomes)
    return plasmids