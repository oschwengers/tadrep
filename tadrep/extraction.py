import logging

import tadrep.config as cfg
import tadrep.io as tio

log = logging.getLogger('EXTRACT')


def extract():
    # get existing json existing_plasmid_dict
    json_output_path = cfg.output_path.joinpath('extraction.json')
    existing_plasmid_dict = tio.load_data(json_output_path)
    
    # update plasmid count
    number_of_plasmids = max([int(id) for id in existing_plasmid_dict.keys()], default=0)
    cfg.verboseprint(f'found {number_of_plasmids} existing plasmids!')
    log.info('found %d existing plasmids in file %s', number_of_plasmids, json_output_path)

    new_plasmids = {}
    number_of_plasmids += 1

    # load sequences
    for input_file in cfg.files_to_extract:

        file_sequences = tio.import_sequences(input_file, sequence=True)
        cfg.verboseprint(f'file: {input_file.name}, sequences: {len(file_sequences)}')
        log.info('file: %s, sequences: %d', input_file.name, len(file_sequences))

        # call genome/draft/plasmid methods
        if(cfg.file_type == 'genome'):
            file_sequences = filter_longest(file_sequences)
        elif(cfg.file_type == 'draft'):
            file_sequences = search_headers(file_sequences)
        
        # add file name and new id to plasmids
        for entry in file_sequences.values():
            entry['file'] = input_file.name
            entry['new_id'] = number_of_plasmids

            new_plasmids[number_of_plasmids] = entry
            number_of_plasmids += 1

    # update existing_plasmid_dict
    existing_plasmid_dict.update(new_plasmids)
    cfg.verboseprint(f'total plasmids found: {len(existing_plasmid_dict)}')
    log.info('total plasmids found: %d', len(existing_plasmid_dict))
    
    # export to json
    tio.export_json(existing_plasmid_dict, json_output_path)


def search_headers(seq_dict):
    filtered_dict = {}
    # search headers for 'plasmid' 'complete', 'circular=true' or custom string
    description_headers = ['plasmid', 'complete', 'circular=true']

    for sequence in seq_dict:
        # test if id contains custom string
        if(cfg.header):
            if(cfg.header in sequence['id']):
                filtered_dict.update(sequence)
                continue
        
        # test description for predefined headers
        if(any(description in sequence['description'] for description in description_headers)):
            filtered_dict.update(sequence)

    return filtered_dict


def filter_longest(seq_dict):
    log.info('discard %d longest sequences from %d entries', cfg.discard, len(seq_dict))
    cfg.verboseprint(f'discard {cfg.discard} longest sequences from {len(seq_dict)} entries')

    # Error if discard too high
    if(cfg.discard > len(seq_dict)):
        log.error('ERROR: can not discard %d sequences from %d present!', cfg.discard, len(seq_dict))
        sys.exit('ERROR: can not discard more sequences than present!')
    
    # sort dict entries by length and discard longest
    filtered_list = sorted(seq_dict.values(), key=lambda x: x['length'], reverse=True)[cfg.discard:]

    # rewrite list into dict
    filtered_dict = {}
    for entry in filtered_list:
        filtered_dict[entry['id']] = entry
    
    return filtered_dict
