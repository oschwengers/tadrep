import logging
import os
import sys
import tempfile
from pathlib import Path


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
    global tmp_path, genome_path, plasmids_path, output_path, prefix, summary_path

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

    for path in args.genome:
        try:
            resolved_path = Path(path).resolve()
            if(not resolved_path.is_file()):
                log.error('genome file %s does not exist!', resolved_path)
                sys.exit(f'ERROR: genome file ({path}) does not exist!')
            if(not os.access(str(resolved_path), os.R_OK)):
                log.error('genome file not readable! path=%s', resolved_path)
                sys.exit(f'ERROR: genome file ({path}) not readable!')
            if(resolved_path.stat().st_size == 0):
                log.error('empty genome file! path=%s', resolved_path)
                sys.exit(f'ERROR: genome file ({path}) is empty!')
            genome_path.append(resolved_path)
        except (OSError, ValueError):
            log.error('provided genome file not valid! path=%s', path)
            sys.exit(f'ERROR: genome file ({path}) not valid!')
        log.info('genome-path=%s', resolved_path)

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
    