import logging
import multiprocessing as mp
import sys
import tempfile
import shutil

from pathlib import Path

import tadrep.utils as tu
import tadrep.io as tio


log = logging.getLogger('CONFIG')

# general setup
# runtime configurations
threads = None
verbose = None
verbose_print = None

# input / output configuration
output_path = None
tmp_path = None
prefix = None

# database setup
force = False
db_type = 'refseq'

# extraction setup
# Input
files_to_extract = None
discard_longest = 1
max_length = None
file_type = None
header = None

# characterize setup
db_local_path = None

# cluster setup
cluster_sequence_identity_threshold = None
cluster_length_threshold= None
skip_cluster = False

# detection setup
# Input
genome_path = []
summary_path = None
db_path = None
db_data = None

# workflow configuration
min_contig_coverage = None
min_contig_identity = None
min_plasmid_coverage = None
min_plasmid_identity = None
gap_sequence_length = None

# multithreading
lock = None
blast_threads = None

# visualize setup
plot_style = 'arrow'
label_color = 'black'
line_width = 0.0
arrow_shaft_ratio = 0.5
size_ratio = 1.0

interval_start = 0.8
interval_number = 10
interval_size = 0.1

label_size = 15
label_rotation = 45
label_hpos = 'center'
label_ha = 'left'

omit_ratio = 1

def setup(args):
    """Test environment and build a runtime configuration."""

    # runtime configurations
    global threads, verbose, verbose_print
    threads = args.threads
    log.info('threads=%i', threads)
    verbose = args.verbose
    log.info('verbose=%s', verbose)
    verbose_print = print if verbose else lambda *a, **k: None

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


def setup_database(args):
    global force, db_type

    force = args.force
    log.info('force = %s', force)
    db_type = args.type
    log.info('db_type = %s', db_type)


def setup_extract(args):
    global files_to_extract, discard_longest, max_length, file_type, header

    if(not args.files):
        log.error('No files provided!')
        sys.exit('ERROR: No input file was provided!')

    files_to_extract = [tu.check_file_permission(file, 'plasmid') for file in args.files]

    discard_longest = args.discard_longest
    if(discard_longest < 0):
        log.error('wrong discard-longest value!')
        sys.exit('ERROR: Wrong parameter value for discard-ongest!')
    
    max_length = args.max_length
    if(max_length < 0):
        log.error('wrong max length value!')
        sys.exit('ERROR: Wrong parameter value for max length!')

    file_type = args.type
    header = args.header
    if(file_type == 'draft' and not header):
        log.debug('No custom header provided!')
        verbose_print('Info: No custom header provided! Only searching for "complete", "circular" and "plasmid"')
    else:
        verbose_print(f'Searching custom header: {header}')
        log.debug('Custom header: %s', header)


def setup_characterize(args):
    global db_local_path

    if(args.database):
        db_global_path = tu.check_file_permission(args.database, 'database')
        db_local_path = output_path.joinpath('db.json')
        shutil.copyfile(db_global_path, db_local_path)
        verbose_print(f'Imported JSON from {db_global_path}')
        log.debug('Copied file from %s to %s', db_global_path, db_local_path)

    if(args.inc_types):
        inc_types_path = tu.check_file_permission(args.inc_types, 'inc-types')
        db_local_path = output_path.joinpath('inc-types.fasta')
        shutil.copyfile(inc_types_path, db_local_path)
        verbose_print(f'Imported inc-types from {inc_types_path}')
        log.debug('Copied file from %s to %s', inc_types_path, db_local_path)


def setup_cluster(args):
    global skip_cluster, cluster_sequence_identity_threshold, cluster_length_threshold

    cluster_sequence_identity_threshold = args.min_sequence_identity / 100
    log.info('cluster-sequence-identity-threshold=%0.3f', cluster_sequence_identity_threshold)

    cluster_length_threshold = args.min_sequence_length_difference
    log.info('cluster-length-threshold=%d', cluster_length_threshold)

    if(args.skip):
        skip_cluster = True
        verbose_print('Skipping clustering')
    log.info('skip_clusters=%s', skip_cluster)


def setup_detect(args):
    # input / output path configurations
    global genome_path, summary_path, db_path, db_data

    if(not args.genome):
        log.error('genome file not provided!')
        sys.exit('ERROR: no genome file was provided!')

    genome_path = [tu.check_file_permission(file, 'genome') for file in args.genome]

    summary_path = output_path.joinpath('summary.tsv')
    log.info('summary_path=%s', summary_path)

    db_path = output_path.joinpath('db.json')
    log.info('db_path=%s', db_path)

    db_data = tio.load_data(db_path)

    if(not db_data):
        log.debug("No data in %s", db_path)
        sys.exit(f"ERROR: No data available in {db_path}")

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


def setup_visualize(args):

    global plot_style, label_color, line_width, arrow_shaft_ratio, size_ratio
    plot_style = args.plot_style
    label_color = args.label_color
    line_width = args.line_width
    arrow_shaft_ratio = args.arrow_shaft_ratio
    size_ratio = args.size_ratio
    log.info('plot_style: %s, label_color: %s, line_width: %f, arrow_shaft_ratio: %f, size_ratio: %f',
    plot_style, label_color, line_width, arrow_shaft_ratio, size_ratio)

    global label_size, label_rotation, label_hpos, label_ha
    label_size = args.label_size
    label_rotation = args.label_rotation
    label_hpos = args.label_hpos
    label_ha = args.label_ha
    log.info('label_size: %d, label_rotation: %d, label_hpos: %s, label_ha: %s',
    label_size, label_rotation, label_hpos, label_ha)

    global interval_start, interval_number, interval_size
    interval_start = args.interval_start / 100
    interval_number = args.interval_number
    interval_size = (1.0 - interval_start) / interval_number
    log.info('Interval: start: %f, count: %d, size: %f', interval_start, interval_number, interval_size)

    global omit_ratio
    omit_ratio = args.omit_ratio
    log.info('omit_ratio: %d', omit_ratio)
