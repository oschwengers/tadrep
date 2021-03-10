import logging
import sys
import shutil
from pathlib import Path

import tadrep
import tadrep.constants as tc
import tadrep.config as cfg
import tadrep.fasta as fasta
import tadrep.utils as tu


def main():
    # parse arguments
    args = tu.parse_arguments()

    ############################################################################
    # Setup logging
    ############################################################################
    try:
        output_path = Path(args.output) if args.output else Path.cwd()
        if(not output_path.exists()):
            output_path.mkdir(parents=True, exist_ok=True)
        output_path = output_path.resolve()
        cfg.output_path = output_path
    except:
        sys.exit(f'ERROR: could not resolve or create output directory ({args.output})!')
    logging.basicConfig(
        filename=str(output_path.joinpath(f'{cfg.prefix}.log')),
        filemode='w',
        format='%(asctime)s - %(levelname)s - %(name)s - %(message)s',
        level=logging.DEBUG if args.verbose else logging.INFO
    )
    log = logging.getLogger('MAIN')
    log.info('version %s', tadrep.__version__)
    log.info('command line: %s', ' '.join(sys.argv))


    ############################################################################
    # Checks and configurations
    # - check parameters and setup global configuration
    # - test database
    # - test binary dependencies
    ############################################################################
    cfg.setup(args)  # check parameters and prepare global configuration
    if(cfg.verbose):
        print(f'TaDReP v{tadrep.__version__}')
        print('Options and arguments:')
        print(f'\tinput: {cfg.genome_path}')
        print(f'\tplasmid(s): {cfg.plasmid_path}')
        print(f'\toutput: {cfg.output_path}')
        print(f'\tprefix: {cfg.prefix}')
        print(f'\ttmp directory: {cfg.tmp_path}')
        print(f'\t# threads: {cfg.threads}')


    ############################################################################
    # Import draft genome contigs
    # - parse contigs in Fasta file
    # - apply contig length filter
    ############################################################################
    print('parse genome sequences...')
    try:
        contigs = fasta.import_contigs(cfg.genome_path)
        log.info('imported sequences=%i', len(contigs))
        print(f'\timported: {len(contigs)}')
    except:
        log.error('wrong genome file format!', exc_info=True)
        sys.exit('ERROR: wrong genome file format!')
    
    # contig length filter (<= default / parameter?)


    ############################################################################
    # Import plasmid sequences
    # - parse contigs in Fasta file
    ############################################################################
    print('parse genome sequences...')
    try:
        contigs = fasta.import_sequences(cfg.genome_path)
        log.info('imported sequences=%i', len(contigs))
        print(f'\timported: {len(contigs)}')
    except:
        log.error('wrong genome file format!', exc_info=True)
        sys.exit('ERROR: wrong genome file format!')


    
    

    ############################################################################
    # Write output files
    # - write comprehensive annotation results as JSON
    # - write optional output files in GFF3/GenBank/EMBL formats
    # - remove temp directory
    ############################################################################
    print('write genome sequences...')
    prefix = f"<prefix>-<plasmid-id>"
    fna_path = cfg.output_path.joinpath(f'{prefix}.fna')
    fasta.export_contigs(hits['contigs'], fna_path, description=True, wrap=True)

    # remove tmp dir
    shutil.rmtree(str(cfg.tmp_path))
    log.debug('removed tmp dir: %s', cfg.tmp_path)


if __name__ == '__main__':
    main()
