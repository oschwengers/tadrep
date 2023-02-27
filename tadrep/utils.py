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
    extraction_parser = subparsers.add_parser('extract', help='Extract unique plasmid sequences')

    arg_group_io = extraction_parser.add_argument_group('Input')
    arg_group_io.add_argument('--type', '-t', action='store', default='genome', choices=['genome', 'plasmid', 'draft'], help='Type of input files')
    arg_group_io.add_argument('--header', action='store', default=None, help='Template for header description inside input files: e.g.: header: ">pl1234" --> --header "pl"')
    arg_group_io.add_argument('--files', '-f', action='store', default=None, nargs="+", help='File path')
    arg_group_io.add_argument('--discard-longest', '-d', action='store', type=int, default=1, help='Discard n longest sequences in output')

    # characterization parser
    characterization_parser = subparsers.add_parser('characterize', help='Identify plasmids with GC content, Inc types, conjugation genes')
    
    arg_group_char = characterization_parser.add_argument_group('Input')
    arg_group_char.add_argument('--json', action='store', default=None, dest='json', help='Import json file from a given database path into working directory')
    arg_group_char.add_argument('--inc-types', action='store', default=None, help='Import inc-types from given path into working directory')

    # clustering parser
    clustering_parser = subparsers.add_parser('cluster', help='Cluster related plasmids')
    
    arg_group_parameters = clustering_parser.add_argument_group('Parameter')
    arg_group_parameters.add_argument('--min-coverage', action='store', type=int, default=90, choices=(1,101), metavar='[1-100]', dest='min_coverage', help='Minimal plasmid coverage (default = 90%%)')
    arg_group_parameters.add_argument('--min-identity', action='store', type=int, default=90, choices=(1,101), metavar='[1-100]', dest='min_coverage', help='Minimal plasmid identity (default = 90%%)')


    # detection parser
    detection_parser = subparsers.add_parser('detect', help='Detect and reconstruct plasmids in draft genomes')

    arg_group_io = detection_parser.add_argument_group('Input / Output')
    arg_group_io.add_argument('--genome', '-g', action='store', default=None, nargs="+", help='Draft genome path')

    arg_group_parameters = detection_parser.add_argument_group('Annotation')
    arg_group_parameters.add_argument('--min-contig-coverage', action='store', type=int, default=90, choices=range(1, 101), metavar='[1-100]', dest='min_contig_coverage', help="Minimal contig coverage (default = 90%%)")
    arg_group_parameters.add_argument('--min-contig-identity', action='store', type=int, default=90, choices=range(1, 101), metavar='[1-100]', dest='min_contig_identity', help="Maximal contig identity (default = 90%%)")
    arg_group_parameters.add_argument('--min-plasmid-coverage', action='store', type=int, default=80, choices=range(1, 101), metavar='[1-100]', dest='min_plasmid_coverage', help="Minimal plasmid coverage (default = 80%%)")
    arg_group_parameters.add_argument('--min-plasmid-identity', action='store', type=int, default=90, choices=range(1, 101), metavar='[1-100]', dest='min_plasmid_identity', help="Minimal plasmid identity (default = 90%%)")
    arg_group_parameters.add_argument('--gap-sequence-length', action='store', type=is_positive, default=10, dest='gap_sequence_length', help="Gap sequence N length (default = 10)")

    # visualization parser
    visualization_parser = subparsers.add_parser('visualize', help='Visualize plasmid coverage of contigs')
    
    arg_group_plot = visualization_parser.add_argument_group('Style')
    arg_group_plot.add_argument('--plotstyle', action='store', default='box', dest='plotstyle', choices=['bigarrow', 'arrow', 'bigbox', 'box', 'bigrbox', 'rbox'], help='Contig representation in plot')
    arg_group_plot.add_argument('--labelcolor', action='store', default='black', dest='labelcolor', help='Contig label color')
    arg_group_plot.add_argument('--linewidth', action='store', default=0.0, type=float, dest='linewidth', help='Contig edge linewidth')
    arg_group_plot.add_argument('--arrow-shaft-ratio', action='store', default=0.5, type=float, dest='arrow_shaft_ratio', help='Size ratio between arrow head and shaft')
    arg_group_plot.add_argument('--size-ratio', action='store', default=1.0, type=float, dest='size_ratio', help='Contig size ratio to track')

    arg_group_label = visualization_parser.add_argument_group('Label')
    arg_group_label.add_argument('--labelsize', action='store', default=15, type=int, dest='labelsize', help='Contig label size')
    arg_group_label.add_argument('--labelrotation', action='store', default=45, type=int, dest='labelrotation', help='Contig label rotation')
    arg_group_label.add_argument('--labelhpos', action='store', default='center', dest='labelhpos', choices=['left', 'center', 'right'], help='Contig label horizontal position')
    arg_group_label.add_argument('--labelha', action='store', default='left', dest='labelha', choices=['left', 'center', 'right'], help='Contig label horizontal alignment')

    arg_group_gradient = visualization_parser.add_argument_group('Gradient')
    arg_group_gradient.add_argument('--interval-start', action='store', type=float, default=80, metavar='[0-100]', dest='interval_start', help='Percentage where gradient should stop')
    arg_group_gradient.add_argument('--number-of-intervals', action='store', type=int, default=10, choices=range(1, 101), metavar='[1-100]', dest='interval_number', help='Number of gradient intervals')

    arg_group_omit = visualization_parser.add_argument_group('Omit')
    arg_group_omit.add_argument('--omit_ratio', action='store', default=1, type=int, choices=range(0, 101), metavar='[0-100]', dest='omit_ratio', help='Omit contigs shorter than X percent of plasmid length from plot')

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
