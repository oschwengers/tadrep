import argparse
import logging
import multiprocessing as mp
import os
import subprocess as sp
import sys

from pathlib import Path

from matplotlib.pyplot import title

import tadrep


log = logging.getLogger('UTILS')

DB_FILES = ['db.tsv', 'db.ndb', 'db.not', 'db.ntf', 'db.nto']

CITATION = 'Schwengers et al. (2022)\nTaDReP: Targeted Detection and Reconstruction of Plasmids.\nGitHub https://github.com/oschwengers/tadrep'


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='TaDReP',
        description='Targeted Detection and Reconstruction of Plasmids',
        epilog=f'Citation:\n{CITATION}\n\nGitHub:\nhttps://github.com/oschwengers/tadrep',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False
    )

    arg_group_general = parser.add_argument_group('General')
    arg_group_general.add_argument('--help', '-h', action='help', help='Show this help message and exit')
    arg_group_general.add_argument('--verbose', '-v', action='store_true', help='Print verbose information')
    arg_group_general.add_argument('--threads', '-t', action='store', type=is_positive, default=mp.cpu_count(), help='Number of threads to use (default = number of available CPUs)')
    arg_group_general.add_argument('--tmp-dir', action='store', default=None, help='Temporary directory to store blast hits')
    arg_group_general.add_argument('--version', action='version', version='%(prog)s ' + tadrep.__version__)

    arg_group_gio = parser.add_argument_group('General Input / Output')
    arg_group_gio.add_argument('--output', '-o', action='store', default=os.getcwd(), help='Output directory (default = current working directory)')
    arg_group_gio.add_argument('--prefix', action='store', default=None, type=not_empty, help='Prefix for all output files (default = None)')

    # add subparser
    subparsers = parser.add_subparsers(dest='subcommand')

    # extraction parser
    extraction_parser = subparsers.add_parser('extraction', help='Extract unique plasmid sequences')

    arg_group_io = extraction_parser.add_argument_group('Input')
    arg_group_io.add_argument('--files', '-f', action='store', default=None, nargs="+", help='File path')
    arg_group_io.add_argument('--drop', '-d', action='store', type=int, default=1, help='Remove n longest circular read(s) in output')

    # characterization parser
    characterization_parser = subparsers.add_parser('characterization', help='Identify plasmids with GC content, Inc types, conjugation genes')

    arg_group_io = characterization_parser.add_argument_group('Input')
    arg_group_io.add_argument('--plasmids', '-p', action='store', default=None, nargs='+', help='Plasmids path')

    # clustering parser
    clustering_parser = subparsers.add_parser('clustering', help='Cluster related plasmids')
    
    arg_group_parameters = clustering_parser.add_argument_group('Parameter')
    arg_group_parameters.add_argument('--min-coverage', action='store', type=int, default=90, choices=(1,101), metavar='[1-100]', dest='min_coverage', help='Minimal plasmid coverage (default = 90%%)')
    arg_group_parameters.add_argument('--min-identity', action='store', type=int, default=90, choices=(1,101), metavar='[1-100]', dest='min_coverage', help='Minimal plasmid identity (default = 90%%)')


    # detection parser
    detection_parser = subparsers.add_parser('detection', help='Detect and reconstruct plasmids in draft genomes')

    arg_group_io = detection_parser.add_argument_group('Input / Output')
    arg_group_io.add_argument('--genome', '-g', action='store', default=None, nargs="+", help='Draft genome path')
    arg_group_io.add_argument('--plasmids', '-p', action='store', default=None, help='Plasmids path')
    arg_group_io.add_argument('--db', action='store', default=None, help='Directory which contains blast database')

    arg_group_parameters = detection_parser.add_argument_group('Annotation')
    arg_group_parameters.add_argument('--min-contig-coverage', action='store', type=int, default=90, choices=range(1, 101), metavar='[1-100]', dest='min_contig_coverage', help="Minimal contig coverage (default = 90%%)")
    arg_group_parameters.add_argument('--min-contig-identity', action='store', type=int, default=90, choices=range(1, 101), metavar='[1-100]', dest='min_contig_identity', help="Maximal contig identity (default = 90%%)")
    arg_group_parameters.add_argument('--min-plasmid-coverage', action='store', type=int, default=80, choices=range(1, 101), metavar='[1-100]', dest='min_plasmid_coverage', help="Minimal plasmid coverage (default = 80%%)")
    arg_group_parameters.add_argument('--min-plasmid-identity', action='store', type=int, default=90, choices=range(1, 101), metavar='[1-100]', dest='min_plasmid_identity', help="Minimal plasmid identity (default = 90%%)")
    arg_group_parameters.add_argument('--gap-sequence-length', action='store', type=is_positive, default=10, dest='gap_sequence_length', help="Gap sequence N length (default = 10)")

    # visualization parser
    visualization_parser = subparsers.add_parser('visualization', help='Visualize plasmid coverage of contigs')
    
    arg_group_plot = visualization_parser.add_argument_group('Plot parameter')
    arg_group_plot.add_argument('--height', action='store', default=0.5,  help='Height per track')
    arg_group_plot.add_argument('--label_y_offset', action='store', default=-3, help='Offset label from track center')
    arg_group_plot.add_argument('--track_edge_distance', action='store', default=1, help='Distance between label and top of previous track')
    arg_group_plot.add_argument('--margin_space', action='store', default=4,  help='Vertical space between plotted elements and plot margin')
    arg_group_plot.add_argument('--label_wrapping', action='store', default=7,  help='Number of characters per contig track')
    arg_group_plot.add_argument('--plasmid_head_modifier', action='store', default=100,  help='Modifier of plasmid head length')
    arg_group_plot.add_argument('--contig_head_modifier', action='store' ,default=200,  help='Modifier of contig head length')

    return parser.parse_args()


def is_positive(value):
    value = int(value)
    if(value < 0):
        raise argparse.ArgumentTypeError(f'{value} is not a positive integer value!')
    return value


def not_empty(string):
    string = str(string).strip()
    if(not string):
        raise argparse.ArgumentTypeError('Prefix cant be an empty string!')
    return string


def run_cmd(cmd_command, tmp_path):
    process = sp.run(
        cmd_command,
        cwd=str(tmp_path),
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )

    if(process.returncode != 0):
        log.debug('command: %s', cmd_command)
        log.debug('stdout=%s, stderr=%s', process.stdout, process.stderr)
        log.warning('command failed! Error-code: %s', process.returncode)
        sys.exit(f'ERROR: {process.stderr}\nError code: {process.returncode}')


def check_file_permission(file, purpose):
    try:
        resolved_path = Path(file).resolve()
        if (not resolved_path.is_file()):
            log.error('%s file %s does not exist!', purpose, resolved_path)
            sys.exit(f'ERROR: {purpose} file ({file}) does not exist!')
        if (not os.access(str(resolved_path), os.R_OK)):
            log.error('%s file not readable! path=%s', purpose, resolved_path)
            sys.exit(f'ERROR: {purpose} file ({file}) not readable!')
        if (resolved_path.stat().st_size == 0):
            log.error('empty %s file! path=%s', purpose, resolved_path)
            sys.exit(f'ERROR: {purpose} file ({file}) is empty!')
    except (OSError, ValueError):
        log.error('provided %s file not valid! path=%s', purpose, file)
        sys.exit(f'ERROR: {purpose} file ({file}) not valid!')
    log.info('%s-path=%s', purpose, resolved_path)
    return resolved_path


def check_db_directory(db_path):
    try:
        resolved_db_path = Path(db_path).resolve()
        if (not resolved_db_path.exists()):
            log.error('Database path does not exist! path=%s', resolved_db_path)
            sys.exit(f'ERROR: database path ({resolved_db_path}) does not exists!')
        for file in DB_FILES:
            file_path = resolved_db_path.joinpath(file)
            if (not file_path.is_file()):
                log.error('Database file missing! file=%s', file_path)
                sys.exit(f'ERROR: Database file missing! File={file_path}')
    except (OSError, ValueError):
        log.error('provided database not valid! path=%s', db_path)
        sys.exit(f'ERROR: provided database {db_path} not valid!')
    return resolved_db_path
