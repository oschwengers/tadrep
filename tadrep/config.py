import logging
import multiprocessing as mp
import sys
import tempfile
import shutil

from pathlib import Path

import tadrep.utils as tu


log = logging.getLogger('CONFIG')

# general setup
# runtime configurations
threads = None
verbose = None
verboseprint = None

# input / output configuration
output_path = None
tmp_path = None
prefix = None

# detection setup
# Input
genome_path = []
summary_path = None

# workflow configuration
min_contig_coverage = None
min_contig_identity = None
min_plasmid_coverage = None
min_plasmid_identity = None
gap_sequence_length = None

# multithreading
lock = None
blast_threads = None

# extraction setup
# Input
files_to_extract = None
discard = 1
file_type = None
header = None

def setup(args):
    """Test environment and build a runtime configuration."""

    # runtime configurations
    global threads, verbose, verboseprint
    threads = args.threads
    log.info('threads=%i', threads)
    verbose = args.verbose
    log.info('verbose=%s', verbose)
    verboseprint = print if verbose else lambda *a, **k: None

    # input / output path configurations
    global tmp_path, output_path, prefix

    if(args.tmp_dir):
        tmp_path = Path(args.tmp_dir)
        if(not tmp_path.exists()):
            log.debug('dedicated temp dir does not exist! tmp-dir=%s', tmp_path)
            sys.exit(f'ERROR: dedicated temp dir ({tmp_path}) does not exist!')
        else:
            log.info('use dedicated temp dir: path=%s', tmp_path)
            tmp_path = Path(tempfile.mkdtemp(dir=str(tmp_path)))
    else:
        tmp_path = Path(tempfile.mkdtemp())
    log.info('tmp-path=%s', tmp_path)

    if(args.prefix):
        prefix = args.prefix
    log.info('output-path=%s', output_path)
    log.info('prefix=%s', prefix)


def setup_detection(args):
    # input / output path configurations
    global genome_path, summary_path

    if(not args.genome):
        log.error('genome file not provided!')
        sys.exit('ERROR: no genome file was provided!')

    genome_path = [tu.check_file_permission(file, 'genome') for file in args.genome]

    summary_path = output_path.joinpath('summary.tsv')
    log.info('summary_path=%s', summary_path)

    # workflow configuration
    global min_contig_coverage, min_contig_identity, min_plasmid_coverage, min_plasmid_identity, gap_sequence_length
    min_contig_coverage = args.min_contig_coverage / 100
    log.info('min-contig-coverage=%0.3f', min_contig_coverage)
    min_contig_identity = args.min_contig_identity / 100
    log.info('min-contig-identity=%0.3f', min_contig_identity)
    min_plasmid_coverage = args.min_plasmid_coverage / 100
    log.info('min-plasmid-coverage=%0.3f', min_plasmid_coverage)
    min_plasmid_identity = args.min_plasmid_identity / 100
    log.info('min-plasmid-identity=%0.3f', min_plasmid_identity)
    gap_sequence_length = args.gap_sequence_length
    log.info('gap-sequence-length=%i', gap_sequence_length)

    # multithreading
    global lock, blast_threads
    lock = mp.Lock()
    blast_threads = threads // len(genome_path)
    if(blast_threads == 0):
        blast_threads = 1
    log.info('blast-threads=%i', blast_threads)


def setup_extraction(args):
    global files_to_extract, discard, file_type, header

    if(not args.files):
        log.error('No files provided!')
        sys.exit('ERROR: No input file was provided!')

    files_to_extract = [tu.check_file_permission(file, 'plasmid') for file in args.files]

    discard = args.discard_longest
    if(discard < 0):
        log.error('Can not drop negative files!')
        sys.exit('ERROR: Can not drop negative files!')

    file_type = args.type
    header = args.header
    if(file_type == 'draft' and not header):
        log.debug('No custom header provided!')
        verboseprint('Info: No custom header provided! Only searching for "complete", "circular" and "plasmid"')
    else:
        verboseprint(f'Searching custom header: {header}')
        log.debug('Custom header: %s', header)


def setup_characterize(args):

    if(args.database):
        json_path = tu.check_file_permission(args.database, 'database')
        target_path = output_path.joinpath('db.json')
        shutil.copyfile(json_path, target_path)
        verboseprint(f'Imported JSON from {json_path}')
        log.debug('Copied file from %s to %s', json_path, target_path)

    if(args.inc_types):
        inc_types_path = tu.check_file_permission(args.inc_types, 'inc-types')
        target_path = output_path.joinpath('inc-types.fasta')
        shutil.copyfile(inc_types_path, target_path)
        verboseprint(f'Imported inc-types from {inc_types_path}')
        log.debug('Copied file from %s to %s', inc_types_path, target_path)
