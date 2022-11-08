import logging

import tadrep.config as cfg
import tadrep.io as tio

log = logging.getLogger('EXTRACT')
verboseprint = print if cfg.verbose else lambda *a, **k: None


def extract():
    # get existing json existing_plasmid_dict
    json_output_path = cfg.output_path.joinpath('extraction.json')
    existing_plasmid_dict = tio.check_existing_JSON(json_output_path)
    
    # update counter    plasmid_count
    plasmid_counter = len(existing_plasmid_dict)
    verboseprint(f'found {plasmid_counter} existing plasmids!')
    log.info('found %d existing plasmids in file %s', plasmid_counter, json_output_path)

    # call genome/draft/plasmid methods
    if(cfg.file_type == 'genome'):
        new_plasmids = genome_extract(plasmid_counter)
    elif(cfg.file_type == 'plasmids'):
        new_plasmids = plasmid_extract(plasmid_counter)
    else:
        new_plasmids = draft_extract(plasmid_counter)

    # update existing_plasmid_dict
    existing_plasmid_dict.update(new_plasmids)
    verboseprint(f'total plasmids found: {len(existing_plasmid_dict)}')
    log.info('total plasmids found: %d', len(existing_plasmid_dict))
    
    # export to json
    tio.export_json(existing_plasmid_dict, json_output_path)


def draft_extract(plasmid_count):
    # remove all sequences without 'circular' 'complete' 'plasmid' or custom header
    log.info('start extraction with "draft" from %d files', len(cfg.files_to_extract))

    new_plasmids = {}

    for input_file in cfg.files_to_extract:

        file_sequences = tio.import_sequences(input_file, sequence=True)
        verboseprint(f'file: {input_file.name}, sequences: {len(file_sequences)}')
        log.info('file: %s, sequences: %d', input_file.name, len(file_sequences))

        file_sequences = search_headers(file_sequences)

        for sequence in file_sequences:
            sequence['file'] = input_file.name

            new_plasmids[plasmid_count] = sequence
            plasmid_count += 1

    return seq_dict


def genome_extract(plasmid_count):
    # read all sequences from input excluding longest
    log.info('start extraction with "genome" from %d files', len(cfg.files_to_extract))
    # find longest X sequences and remove them
    new_plasmids = {}

    for input_file in cfg.files_to_extract:

        file_sequences = tio.import_sequences(input_file, sequence=True)
        verboseprint(f'file: {input_file.name}, sequences: {len(file_sequences)}')
        log.info('file: %s, sequences: %d', input_file.name, len(file_sequences))

        file_sequences = filter_longest(file_sequences)

        for sequence in file_sequences:
            sequence['file'] = input_file.name

            new_plasmids[plasmid_count] = sequence
            plasmid_count += 1

    return new_plasmids


def plasmid_extract(plasmid_count):
    # read all sequences from input file without filtering
    log.info('start extraction with "plasmid" from %d files', len(cfg.files_to_extract))
    # create new empty dict
    new_plasmids = {}

    # for each file
    for input_file in cfg.files_to_extract:
        # read file     tio.import_sequences sequences=true
        file_plasmids = tio.import_sequences(input_file, sequences=True)
        verboseprint(f'file: {input_file.name}, sequences: {len(file_plasmids)}')
        log.info('file: %s, sequences: %d', input_file.name, len(file_plasmids))
        
        for plasmid in file_plasmids:
            # add filename to each sequence in dict
            plasmid['file'] = input_file.name

            # add each plasmid to new dict with plasmid_count as key
            new_plasmids[plasmid_count] = plasmid
            plasmid_count += 1
        
    return new_plasmids


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
    log.info('drop %d longest sequences from %d entries', cfg.drop, len(seq_dict))
    seq_dict = sorted(seq_dict, key=lambda x: x['length'])[cfg.drop:]
    return seq_dict
