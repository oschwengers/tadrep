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
verboseprint = None

# input / output configuration
output_path = None
tmp_path = None
prefix = None

# detection setup
# Input
genome_path = []
summary_path = None
db_path = None
db_data = None
blastdb_path = None

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

# visualize setup
plotstyle = 'arrow'
labelcolor = 'black'
facecolor = 'orange'
linewidth = 0.0
arrow_shaft_ratio = 0.5
size_ratio = 1.0

labelsize = 15
labelrotation = 45
labelhpos = 'center'
labelha = 'left'

omit_ratio = 1

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
    global genome_path, summary_path, db_path, db_data, blastdb_path

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

    blastdb_path = db_data.get('db_path', '')
    log.debug("BlastDB path = %s", blastdb_path)

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

    if(args.json):
        json_path = tu.check_file_permission(args.json, 'database')
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


def setup_visualize(args):

    global plotstyle, labelcolor, facecolor, linewidth, arrow_shaft_ratio, size_ratio
    plotstyle = args.plotstyle
    labelcolor = args.labelcolor
    facecolor = args.facecolor
    linewidth = args.linewidth
    arrow_shaft_ratio = args.arrow_shaft_ratio
    size_ratio = args.size_ratio
    log.info('plotstyle: %s, labelcolor: %s, facecolor: %s, linewdith: %f, arrow_shaft_ratio: %f, size_ratio: %f',
    plotstyle, labelcolor, facecolor, linewidth, arrow_shaft_ratio, size_ratio)

    global labelsize, labelrotation, labelhpos, labelha
    labelsize = args.labelsize
    labelrotation = args.labelrotation
    labelhpos = args.labelhpos
    labelha = args.labelha
    log.info('labelsize: %d, labelrotation: %d, labelhpos: %s, labelha: %s',
    labelsize, labelrotation, labelhpos, labelha)

    global omit_ratio
    omit_ratio = args.omit_ratio
    log.info('omit_ratio: %d', omit_ratio)

