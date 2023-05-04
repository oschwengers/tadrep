import argparse
import logging
import os

import tadrep.utils as tu


log = logging.getLogger('UTILS')


REFSEQ = 'refseq'
PLSDB = 'plsdb'

def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='TaDReP_DB',
        description='Download and create database for TaDReP',
        add_help=False
    )

    arg_group_io = parser.add_argument_group('Input / Output')
    arg_group_io.add_argument('--type', action='store', default='refseq', choices=[REFSEQ, PLSDB], type=str.lower, help="External DB to import (default = 'refseq')")
    arg_group_io.add_argument('--output', '-o', action='store', default=os.getcwd(), help='Output directory for database files (default = current working directory)')

    arg_group_general = parser.add_argument_group('General')
    arg_group_general.add_argument('--help', '-h', action='help', help='Show this help message and exit')
    arg_group_general.add_argument('--verbose', '-v', action='store_true', help='Print verbose information')
    arg_group_general.add_argument('--force', '-f', action='store_true', help='Force download and new setup of database')
    return parser.parse_args()


def reverse_database(fasta_path, db_path, tmp_path):
    log.debug('reverse database: input-path=%s, output-path=%s', db_path, fasta_path)
    cmd_reverse = [
        'blastdbcmd',
        '-entry', 'all',
        '-db', str(db_path),
        '-out', str(fasta_path)
    ]
    log.debug('cmd=%s', cmd_reverse)
    tu.run_cmd(cmd_reverse, tmp_path)
