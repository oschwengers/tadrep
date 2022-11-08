import logging

import tadrep.config as cfg
import tadrep.io as tio

log = logging.getLogger('EXTRACT')


def extract():
    # get existing json existing_plasmid_dict
    json_output_path = cfg.output_path.joinpath('extraction.json')
    existing_plasmid_dict = tio.check_existing_JSON(json_output_path)
    
    # update counter    plasmid_count
    plasmid_counter = len(existing_plasmid_dict)

    # call genome/draft/plasmid methods
    if(cfg.file_type == 'genome'):
        new_plasmids = genome_extract(plasmid_counter)
    elif(cfg.file_type == 'plasmids'):
        new_plasmids = plasmids_extract(plasmid_counter)
    else:
        new_plasmids = draft_extract(plasmid_counter)

    # update existing_plasmid_dict
    existing_plasmid_dict.update(new_plasmids)
    
    # export to json
    tio.export_json(existing_plasmid_dict, json_output_path)


def draft_extract(plasmid_count):
    # remove all sequences without 'circular' 'complete' 'plasmid' or custom header
    return seq_dict


def genome_extract(seq_dict):
    # read all sequences from input excluding longest
    # find longest X sequences and remove them
    return seq_dict


def plasmids_extract(seq_dict):
    # read all sequences from input file without filtering
    return seq_dict


def search_headers(seq_dict):
    # search headers for 'plasmid' 'complete' or custom string
    return seq_dict


def filter_longest(seq_dict):
    # filter out longest sequences
    return seq_dict
