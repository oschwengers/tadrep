import logging
import sys

import tadrep.io as tio
import tadrep.config as cfg

log = logging.getLogger('CHARACTERIZE')


def characterize():
    # load json from output path
    json_path = cfg.output_path.joinpath('extraction.json')
    existing_data = tio.load_data(json_path)
    plasmids = existing_data.get('extraction', {})

    # check if data is available
    if(not existing_data or not plasmids):
        log.error('Failed to load data from %s', json_path)
        sys.exit(f'ERROR: Failed to load data from {json_path}!')

    for plasmid in plasmids:
        # set length
        plasmid['length'] = len(plasmid['sequence'])

        # set gc_content
        plasmid['gc_content'] = gc_content(plasmid['sequence'])
        # calculate INC_types (Platon)

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


def inc_types():
    pass
