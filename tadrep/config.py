import logging
import os
import sys
import tempfile
from pathlib import Path

import tadrep.utils as tu


log = logging.getLogger('CONFIG')


# runtime configurations
threads = None
verbose = None

# input / output configuration
genome_path = []
plasmids_path = None
output_path = None
tmp_path = None
summary_path = None
prefix = None
database_path = None

# workflow configuration
min_contig_coverage = None
min_contig_identity = None
min_plasmid_coverage = None
min_plasmid_identity = None
gap_sequence_length = None


def setup(args):
    """Test environment and build a runtime configuration."""
    # runtime configurations
    global threads, verbose
    threads = args.threads
    log.info('threads=%i', threads)
    verbose = args.verbose
    log.info('verbose=%s', verbose)

    # input / output path configurations
    global tmp_path, genome_path, plasmids_path, output_path, prefix, summary_path, database_path

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

    if(not args.genome):
        log.error('genome file not provided!')
        sys.exit('ERROR: no genome file was provided!')

    genome_path = [tu.check_file_permission(file, 'genome') for file in args.genome]

    if(args.plasmids):
        plasmids_path = tu.check_file_permission(args.plasmids, 'plasmids')
    elif(args.database):
        database_path = tu.check_db_directory(args.database)
        log.info('database path=%s', database_path)
    else:
        log.error('no plasmid file or database was provided!')
        sys.exit('ERROR: neither plasmid file nor database was provided!')

    if(args.prefix):
        prefix = args.prefix
    log.info('plasmids-path=%s', plasmids_path)
    log.info('output-path=%s', output_path)
    log.info('prefix=%s', prefix)

    # workflow configuration
    global min_contig_coverage, min_contig_identity, min_plasmid_coverage, min_plasmid_identity, gap_sequence_length
    min_contig_coverage = args.min_contig_coverage / 100
    log.info('min-contig-coverage=%s', min_contig_coverage)
    min_contig_identity = args.min_contig_identity / 100
    log.info('min-contig-identity=%s', min_contig_identity)
    min_plasmid_coverage = args.min_plasmid_coverage / 100
    log.info('min-plasmid-coverage=%s', min_plasmid_coverage)
    min_plasmid_identity = args.min_plasmid_identity / 100
    log.info('min-plasmid-identity=%s', min_plasmid_identity)
    gap_sequence_length = args.gap_sequence_length
    log.info('gap-sequence-length=%s', gap_sequence_length)
