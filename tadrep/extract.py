import logging
import sys

import tadrep.config as cfg
import tadrep.io as tio

log = logging.getLogger('EXTRACT')


def extract():
    # get existing json existing_plasmid_dict
    json_output_path = cfg.output_path.joinpath('db.json')
    db_data = tio.load_data(json_output_path)
    plasmid_dict = db_data.get('plasmids', {})      # are previous sequences available
    file_list = db_data.get('files', [])              # which files were already extracted from
    
    # update plasmid count
    number_of_plasmids = len(plasmid_dict.keys())
    cfg.verbose_print(f'Loaded {number_of_plasmids} previously extracted plasmids')
    log.info('Loaded %d previously extracted plasmids from file %s', number_of_plasmids, json_output_path)

    new_plasmids = {}

    # load sequences
    for input_file in cfg.files_to_extract:
        # check if file was already extracted
        if(str(input_file) in file_list):
            cfg.verbose_print(f'Skipping {input_file.name}, already extracted from!')
            log.info('Skipping %s, already extracted from!', input_file)
            continue
    
        # add complete filepath to extracted files
        file_list.append(str(input_file))

        # import sequence
        imported_sequences = tio.import_sequences(input_file, sequence=True)
        imported_sequences = imported_sequences.values()
        log.info('File: %s, sequences: %d', input_file.name, len(imported_sequences))

        # call genome/draft/plasmid methods
        extracted_plasmids = []
        if(cfg.file_type == 'genome'):
            extracted_plasmids = filter_longest(imported_sequences)
        elif(cfg.file_type == 'draft'):
            extracted_plasmids = filter_by_header(imported_sequences)
            extracted_plasmids = [plasmid for plasmid in extracted_plasmids if plasmid['length'] <= cfg.max_length]  # apply length filter
        else:  # plasmid
            extracted_plasmids = imported_sequences
        
        # add file name and new id to plasmids
        cfg.verbose_print(f'File: {input_file.name}, detected plasmids: {len(extracted_plasmids)}')
        log.info('File: %s, plasmids detected: %d', input_file.name, len(extracted_plasmids))
        for plasmid in extracted_plasmids:
            plasmid['file'] = input_file.name
            new_plasmids[plasmid['id']] = plasmid

    # update existing_plasmid_dict
    plasmid_dict.update(new_plasmids)
    cfg.verbose_print(f'New plasmids extracted: {len(new_plasmids)}')
    cfg.verbose_print(f'Total plasmids extracted: {len(plasmid_dict)}')
    log.info('Total plasmids extracted: %d', len(plasmid_dict))
    
    # export to json
    db_data['plasmids'] = plasmid_dict
    db_data['files'] = file_list
    tio.export_json(db_data, json_output_path)


def filter_by_header(sequences):
    # search headers for 'plasmid' 'complete', 'circular=true' or custom string
    description_headers = ['plasmid', 'complete', 'circular=true']

    filtered_sequences = []
    for sequence in sequences:
        # test if id contains custom string
        if(cfg.header):
            if(cfg.header in sequence['id']):
                filtered_sequences.append(sequence)
        else:
            # test description for predefined headers
            if(any(description in sequence['description'].lower() for description in description_headers)):
                filtered_sequences.append(sequence)

    return filtered_sequences


def filter_longest(plasmids):
    log.info('Discard %d longest sequences from %d entries', cfg.discard_longest, len(plasmids))
    cfg.verbose_print(f'Discard {cfg.discard_longest} longest sequences from {len(plasmids)} entries')

    # Error if discard too high
    if(cfg.discard_longest > len(plasmids)):
        log.error('Can not discard %d sequences from %d present!', cfg.discard_longest, len(plasmids))
        sys.exit('ERROR: Can not discard more sequences than present!')
    
    # sort dict entries by length and discard longest
    filtered_plasmids = sorted(plasmids, key=lambda x: x['length'], reverse=True)[cfg.discard_longest:]

    return filtered_plasmids
