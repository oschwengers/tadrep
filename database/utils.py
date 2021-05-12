import logging
import argparse
import os

import tadrep.fasta as tf
import tadrep.utils as tu

log = logging.getLogger('UTILS')


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='TaDReP - Database',
        description='Download and create database for TaDReP',
        add_help=False
    )

    arg_group_io = parser.add_argument_group('Input / Output')
    arg_group_io.add_argument('--db', action='store', default='RefSeq', choices=['RefSeq', 'PLSDB'], help="Extern DB to import (default = 'RefSeq')")
    arg_group_io.add_argument('--output', '-o', action='store', default=os.getcwd(), help='Output directory for database files (default = current working directory)')

    arg_group_general = parser.add_argument_group('General')
    arg_group_general.add_argument('--help', '-h', action='help', help='Show this help message and exit')
    arg_group_general.add_argument('--verbose', '-v', action='store_true', help='Print verbose information')
    arg_group_general.add_argument('--force', '-f', action='store_true', help='Force download and new setup of database')
    return parser.parse_args()


def create_tsv(fasta_path, tsv_output_path):
    log.info('Create .tsv: input=%s, output=%s', fasta_path, tsv_output_path)
    plasmids = tf.import_sequences(fasta_path)
    with tsv_output_path.open('w') as fh_out:
        for plasmid_id, plasmid in plasmids.items():
            fh_out.write(f"{plasmid['id']}\t{plasmid['description']}\t{plasmid['length']}\n")


def create_blast_db(output_path, fasta_tmp_path, tmp_path):

    cmd_blastdb = [
        'makeblastdb',
        '-in', str(fasta_tmp_path),
        '-dbtype', 'nucl',
        '-out', str(output_path),
        '-parse_seqids'
    ]
    log.debug('cmd_blastdb=%s', cmd_blastdb)
    tu.run_cmd(cmd_blastdb, tmp_path)
