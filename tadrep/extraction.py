import logging

import tadrep.config as cfg
import tadrep.io as tio

log = logging.getLogger('EXTRACT')
verboseprint = print if cfg.verbose else lambda *a, **k: None


def extract():
    # get existing json existing_plasmid_dict
    json_output_path = cfg.output_path.joinpath('extraction.json')
    existing_plasmid_dict = tio.load_data(json_output_path)
    
    # update plasmid count
    number_of_plasmids = max([id for id in existing_plasmid_dict.keys()])
    verboseprint(f'found {number_of_plasmids} existing plasmids!')
    log.info('found %d existing plasmids in file %s', number_of_plasmids, json_output_path)

    new_plasmids = {}

    # load sequences
    for input_file in cfg.files_to_extract:

        file_sequences = tio.import_sequences(input_file, sequences=True)
        verboseprint(f'file: {input_file.name}, sequences: {len(file_sequences)}')
        log.info('file: %s, sequences: %d', input_file.name, len(file_sequences))

        # call genome/draft/plasmid methods
        if(cfg.file_type == 'genome'):
            file_sequences = filter_longest(file_sequences)
        elif(cfg.file_type == 'draft'):
            file_sequences = search_headers(file_sequences)
        
        # add file name and new id to plasmids
        for sequence in file_sequences:
            sequence['file'] = input_file.name
            sequence['new_id'] = plasmid_count

            new_plasmids[plasmid_count] = sequence
            plasmid_count += 1

    # update existing_plasmid_dict
    existing_plasmid_dict.update(new_plasmids)
    verboseprint(f'total plasmids found: {len(existing_plasmid_dict)}')
    log.info('total plasmids found: %d', len(existing_plasmid_dict))
    
    # export to json
    tio.export_json(existing_plasmid_dict, json_output_path)


def search_headers(seq_dict):
    filtered_dict = {}
    # search headers for 'plasmid' 'complete' or custom string
    standard_headers = ['plasmid', 'complete', 'circular']
    if(cfg.header):
        standard_headers.append(cfg.header)

    log.info('searching %s', standard_headers)
    
    for entry in seq_dict:
        if(any(substring in entry['description'] for substring in standard_headers)):
            filtered_dict.update(entry)
    
    log.info('found %d matching sequences', len(filtered_dict))
    return filtered_dict


def filter_longest(seq_dict):
    log.info('drop %d longest sequences from %d entries', cfg.discard, len(seq_dict))
    seq_dict = sorted(seq_dict, key=lambda x: x['length'])[cfg.discard:]
    return seq_dict
