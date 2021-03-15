
import logging
import os
import sys
import tempfile
from pathlib import Path

log = logging.getLogger('CONFIG')

# runtime configurations
env = os.environ.copy()
threads = None
verbose = None

# input / output configuration
genome_path = None
plasmids_path = None
output_path = None
tmp_path = None

# workflow configuration
min_contig_coverage = 0.9
min_contig_identity = 0.9
min_plasmid_coverage = 0.8
min_plasmid_identity = 0.8
gap_sequence_length = 10


def setup(args):
    """Test environment and build a runtime configuration."""
    # runtime configurations
    global env, threads, verbose
    threads = args.threads
    log.info('threads=%i', threads)
    verbose = args.verbose
    log.info('verbose=%s', verbose)

    # input / output path configurations
    global tmp_path, genome_path, plasmids_path, output_path

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

    try:
        genome_path = Path(args.genome).resolve()
        if(not os.access(str(genome_path), os.R_OK)):
            log.error('genome file not readable! path=%s', genome_path)
            sys.exit(f'ERROR: genome file ({genome_path}) not readable!')
        if(genome_path.stat().st_size == 0):
            log.error('empty genome file! path=%s', genome_path)
            sys.exit(f'ERROR: genome file ({genome_path}) is empty!')
    except:
        log.error('provided genome file not valid! path=%s', args.genome)
        sys.exit(f'ERROR: genome file ({args.genome}) not valid!')
    log.info('genome-path=%s', genome_path)

    try:
        plasmids_path = Path(args.plasmids).resolve()
        if(not os.access(str(plasmids_path), os.R_OK)):
            log.error('plasmids file not readable! path=%s', plasmids_path)
            sys.exit(f'ERROR: plasmids file ({plasmids_path}) not readable!')
        if(plasmids_path.stat().st_size == 0):
            log.error('empty plasmids file! path=%s', plasmids_path)
            sys.exit(f'ERROR: plasmids file ({plasmids_path}) is empty!')
    except:
        log.error('provided plasmids file not valid! path=%s', args.plasmids)
        sys.exit(f'ERROR: plasmids file ({args.plasmids}) not valid!')
    log.info('plasmids-path=%s', plasmids_path)
    log.info('output-path=%s', output_path)

    # workflow configuration
    global min_contig_coverage, min_contig_identity, min_plasmid_coverage, min_plasmid_identity, gap_sequence_length
    if(args.min_contig_coverage):
        min_contig_coverage = int(args.min_contig_coverage) / 100
    log.info('min-contig-coverage=%s', min_contig_coverage)
    if(args.min_contig_identity):
        min_contig_identity = int(args.min_contig_identity) / 100
    log.info('min-contig-identity=%s', min_contig_identity)
    if(args.min_plasmid_coverage):
        min_plasmid_coverage = int(args.min_plasmid_coverage) / 100
    log.info('min-plasmid-coverage=%s', min_plasmid_coverage)
    if(args.min_plasmid_identity):
        min_plasmid_identity = int(args.min_plasmid_identity) / 100
    log.info('min-plasmid-identity=%s', min_plasmid_identity)
    gap_sequence_length = args.gap_sequence_length
    log.info('gap-sequence-length=%s', gap_sequence_length)
    