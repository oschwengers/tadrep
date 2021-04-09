import argparse
import logging
import multiprocessing as mp
import os

import tadrep
import tadrep.constants as bc
import tadrep.config as cfg

log = logging.getLogger('UTILS')


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='TaDReP',
        description='Targeted Detection and Reconstruction of Plasmids',
        add_help=False
    )

    arg_group_io = parser.add_argument_group('Input / Output')
    arg_group_io.add_argument('--genome', '-g', action='store', default=None, nargs="+", help='Draft genome path (default = <tadrep_path>/db)')
    arg_group_io.add_argument('--plasmids', '-p', action='store', default=None, help='Plasmids path (default = <tadrep_path>/db)')
    arg_group_io.add_argument('--output', '-o', action='store', default=os.getcwd(), help='Output directory (default = current working directory)')
    arg_group_io.add_argument('--tmp-dir', action='store', default=None, help='Temporary directory to store temporary files with blast hits')
    arg_group_io.add_argument('--prefix', action='store', default=None, help='Prefix for all output files')
    
    arg_group_parameters = parser.add_argument_group('Annotation')
    arg_group_parameters.add_argument('--min-contig-coverage', action='store', type=int, default=40, dest='min_contig_coverage', help="Minimal contig coverage (default = 40%%)")
    arg_group_parameters.add_argument('--min-contig-identity', action='store', type=int, default=90, dest='min_contig_identity', help="Maximal contig identity (default = 90%%)")
    arg_group_parameters.add_argument('--min-plasmid-coverage', action='store', type=int, default=80, dest='min_plasmid_coverage', help="Minimal plasmid coverage (default = 80%%)")
    arg_group_parameters.add_argument('--min-plasmid-identity', action='store', type=int, default=90, dest='min_plasmid_identity', help="Minimal plasmid identity (default = 90%%)")
    arg_group_parameters.add_argument('--gap-sequence-length', action='store', type=int, default=10, dest='gap_sequence_length', help="Gap sequence N length (default = 10)")

    arg_group_general = parser.add_argument_group('General')
    arg_group_general.add_argument('--help', '-h', action='help', help='Show this help message and exit')
    arg_group_general.add_argument('--verbose', '-v', action='store_true', help='Print verbose information')
    arg_group_general.add_argument('--threads', '-t', action='store', type=int, default=mp.cpu_count(), help='Number of threads to use (default = number of available CPUs)')
    arg_group_general.add_argument('--version', action='version', version='%(prog)s ' + tadrep.__version__)
    return parser.parse_args()