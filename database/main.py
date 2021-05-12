import logging
import sys
import shutil
from pathlib import Path
import tempfile

import tadrep
import database.utils as du
import database.refseq as dr
import database.plsdb as dp


def main():
    args = du.parse_arguments()
    try:
        output_path = Path(args.output) if args.output else Path.cwd()
        if(not output_path.exists()):
            output_path.mkdir(parents=True, exist_ok=True)
        output_path = output_path.resolve()
    except:
        sys.exit(f'ERROR: could not resolve output directory!')

    logging.basicConfig(
        filename=str(output_path.joinpath('tadrepdb.log')),
        filemode='w',
        format='%(asctime)s - %(levelname)s - %(name)s - %(message)s',
        level=logging.DEBUG if args.verbose else logging.INFO
    )
    log = logging.getLogger('MAIN')
    log.info('version %s', tadrep.__version__)
    log.info('command line: %s', ' '.join(sys.argv))
    log.info('output_path=%s', output_path)

    verbose = args.verbose
    log.info('verbose=%s', verbose)
    force = args.force
    log.info('force=%s', force)

    selected_db = args.db
    log.info('Selected database: %s', selected_db)

    tmp_path = Path(tempfile.mkdtemp())
    log.info('tmp-path=%s', tmp_path)

    if(verbose):
        print(f'TaDReP Database Creation v{tadrep.__version__}')
        print('Options and arguments:')
        print(f'\tdatabase: {selected_db}')
        print(f"\toutput: {output_path}")
        print(f'\ttmp directory: {tmp_path}')

    if(selected_db == 'RefSeq'):
        db_name = 'refseq'
    else:
        db_name = 'plsdb'

    db_output_path = output_path.joinpath(f'{db_name}')
    if(db_output_path.exists()):
        if(force):
            shutil.rmtree(db_output_path)
            log.info('Directory removed: path=%s', db_output_path)
        else:
            shutil.rmtree(str(tmp_path))
            log.debug('removed tmp dir: %s', tmp_path)
            log.debug('Database directory already exists! path=%s', db_output_path)
            sys.exit(f'ERROR: Directory "{db_output_path.stem}" already exists in {output_path} !')

    db_output_path.mkdir(parents=True, exist_ok=True)
    log.info('Directory created: path=%s', db_output_path)

    fasta_tmp_path = tmp_path.joinpath(f'{db_name}.fna')
    log.info('Merged fasta file: path=%s', fasta_tmp_path)

    print('Download starting...')
    if(selected_db == 'RefSeq'):
        dr.download_database(fasta_tmp_path)
    else:
        dp.download_database(tmp_path)
        print('Fetch database information...')
        dp.reverse_database(fasta_tmp_path, tmp_path)
    print('Create blast database...')
    db_file_name = db_output_path.joinpath('db')
    log.info('Database files: name=%s', db_file_name)
    du.create_blast_db(db_file_name, fasta_tmp_path, tmp_path)
    
    
    tsv_output_path = db_output_path.joinpath('db.tsv')
    log.info('TSV file: name=%s, path=%s', tsv_output_path.stem, tsv_output_path)
    du.create_tsv(fasta_tmp_path, tsv_output_path)
    
    print(f'Database successfully created\nDatabase path: {db_output_path}')

    shutil.rmtree(str(tmp_path))
    log.debug('removed tmp dir: %s', tmp_path)


if(__name__ == "__main__"):
    main()
