import argparse
import logging
import multiprocessing as mp
import os
import sys
import re
import subprocess as sp
from pathlib import Path

import tadrep


log = logging.getLogger('UTILS')
FILES_V4 = ['db.nhr', 'db.nin', 'db.nsq', 'db.tsv']
FILES_V5 = ['db.not', 'db.ntf', 'db.nto', 'db.ndb']

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
    arg_group_io.add_argument('--database', '-db', action='store', default=None, help='Path to directory which contains database files')
    
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


def run_cmd(cmd_command, tmp_path):
    process = sp.run(
        cmd_command,
        cwd=str(tmp_path),
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )

    if(process.returncode != 0):
        log.debug('blast command: %s', cmd_command)
        log.debug('stdout=%s, stderr=%s', process.stdout, process.stderr)
        log.warning('Blast command failed! Error-code: %s', process.returncode)
        sys.exit(f'ERROR: {process.stderr}\nBlast command error! Error code: {process.returncode}')


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
        files = FILES_V4 + FILES_V5 if check_db_version(resolved_db_path) else FILES_V4
        for file in files:
            file_path = resolved_db_path.joinpath(file)
            if (not file_path.is_file()):
                log.error('Database file missing! file=%s', file_path)
                sys.exit(f'ERROR: Database file missing! File={file_path}')
    except (OSError, ValueError):
        log.error('provided database not valid! path=%s', db_path)
        sys.exit(f'ERROR: provided database {db_path} not valid!')
    return resolved_db_path


def check_db_version(db_path):
    cmd_version = [
        'blastdbcmd',
        '-db', str(db_path.joinpath('db')),
        '-info'
    ]
    log.info('Blastdb cmd: %s', cmd_version)
    regex = r'(Version:\s)(\d+)'
    try:
        db_output = sp.check_output(cmd_version, stderr=sp.STDOUT, universal_newlines=True)
        log.debug('Blastdb cmd output: %s', db_output)
        result = re.search(regex, str(db_output))
        log.info('Regex result: %s', result)
    except FileNotFoundError:
        log.error('Dependency not found! tool=%s', cmd_version[0])
        sys.exit(f'ERROR: Command {cmd_version[0]} not found or not executable!')
    except sp.CalledProcessError as e:
        log.error('Dependency failed! tool=%s', cmd_version[0])
        log.warning('Subprocess error: %s', e.stdout)
        sys.exit(f'ERROR: {cmd_version[0]} could not be executed!\n{e.stdout}')

    if(not result):
        return False
    log.info('Database version: %s', result.group(2))
    if(int(result.group(2)) == 5):
        return True
    return False
