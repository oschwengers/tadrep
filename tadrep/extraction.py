import logging

import tadrep.config as cfg
import tadrep.io as tio

log = logging.getLogger('EXTRACT')


def extract():
    # get existing json existing_plasmid_dict
    # update counter    plasmid_count

    # load files
        # call genome/draft/plasmid methods
        # update existing_plasmid_dict
    # export to json
    pass


def draft_extract(seq_dict):
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
