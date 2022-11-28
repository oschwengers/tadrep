import logging
import sys

import tadrep.config as cfg
import tadrep.io as tio

log = logging.getLogger('EXTRACT')


def extract():
    # get existing json existing_plasmid_dict
    json_output_path = cfg.output_path.joinpath('db.json')
    existing_data = tio.load_data(json_output_path)
    plasmid_dict = existing_data.get('plasmids', {})      # are previous sequences available
    file_list = existing_data.get('files', [])              # which files were already extracted from
    
    # update plasmid count
    number_of_plasmids = max([int(id) for id in plasmid_dict.keys()], default=0)
    cfg.verboseprint(f'Loaded {number_of_plasmids} previously extracted plasmids!')
    log.info('Loaded %d previously extracted plasmids from file %s', number_of_plasmids, json_output_path)

    new_plasmids = {}
    number_of_plasmids += 1

    # load sequences
    for input_file in cfg.files_to_extract:
        # check if file was already extracted
        if(str(input_file) in file_list):
            cfg.verboseprint(f'Skipping {input_file.name}, already extracted from!')
            log.info('Skipping %s, already extracted from!', input_file)
            continue
    
        # add complete filepath to extracted files
        file_list.append(str(input_file))

        # import sequence
        plasmids = tio.import_sequences(input_file, sequence=True)
        cfg.verboseprint(f'File: {input_file.name}, sequences: {len(plasmids)}')
        log.info('File: %s, sequences: %d', input_file.name, len(plasmids))

        # call genome/draft/plasmid methods
        if(cfg.file_type == 'genome'):
            plasmids = filter_longest(plasmids)
        elif(cfg.file_type == 'draft'):
            plasmids = search_headers(plasmids)
        
        # add file name and new id to plasmids
        cfg.verboseprint(f'File: {input_file.name}, detected plasmids: {len(plasmids)}')
        log.info('File: %s, plasmids detected: %d', input_file.name, len(plasmids))
        for plasmid in plasmids.values():
            plasmid['file'] = input_file.name
            plasmid['new_id'] = number_of_plasmids

            new_plasmids[number_of_plasmids] = plasmid
            number_of_plasmids += 1

    # update existing_plasmid_dict
    plasmid_dict.update(new_plasmids)
    cfg.verboseprint(f'Total plasmids detected: {len(plasmid_dict)}')
    log.info('Total plasmids detected: %d', len(plasmid_dict))
    
    # export to json
    existing_data['plasmids'] = plasmid_dict
    existing_data['files'] = file_list
    tio.export_json(existing_data, json_output_path)


def search_headers(seq_dict):
    filtered_dict = {}
    # search headers for 'plasmid' 'complete', 'circular=true' or custom string
    description_headers = ['plasmid', 'complete', 'circular=true']

    for seq_id, sequence in seq_dict.items():
        # test if id contains custom string
        if(cfg.header):
            if(cfg.header in sequence['id']):
                filtered_dict[seq_id] = sequence
                continue
        else:
            # test description for predefined headers
            if(any(description in sequence['description'] for description in description_headers)):
                filtered_dict[seq_id] = sequence

    return filtered_dict


def filter_longest(seq_dict):
    log.info('Discard %d longest sequences from %d entries', cfg.discard, len(seq_dict))
    cfg.verboseprint(f'Discard {cfg.discard} longest sequences from {len(seq_dict)} entries')

    # Error if discard too high
    if(cfg.discard > len(seq_dict)):
        log.error('Can not discard %d sequences from %d present!', cfg.discard, len(seq_dict))
        sys.exit('ERROR: Can not discard more sequences than present!')
    
    # sort dict entries by length and discard longest
    filtered_list = sorted(seq_dict.values(), key=lambda x: x['length'], reverse=True)[cfg.discard:]

    # rewrite list into dict
    filtered_dict = {}
    for entry in filtered_list:
        filtered_dict[entry['id']] = entry
    
    return filtered_dict
