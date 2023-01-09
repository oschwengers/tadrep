import logging
import shutil
import sys
import tempfile

from pathlib import Path

import tadrep
import tadrep.utils as tu
import tadrep.io as tio
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
    log.info('output-path=%s', output_path)

    verbose = args.verbose
    log.info('verbose=%s', verbose)
    force = args.force
    log.info('force=%s', force)

    selected_db = args.type
    log.info('selected database: %s', selected_db)
    input_files = args.files
    log.info('provided files: %s', input_files)

    update_path = args.database
    log.info('update database: %s', update_path)

    tmp_path = Path(tempfile.mkdtemp())
    log.info('tmp-path=%s', tmp_path)

    if(not update_path):
        if(selected_db == du.REFSEQ):
            db_name = du.REFSEQ
        elif(selected_db == du.PLSDB):
            db_name = du.PLSDB
        else:
            db_name = du.CUSTOM
    else:
        update_path = tu.check_db_directory(update_path)
        output_path = update_path.parent
        db_name = update_path.name

    if(verbose):
        print(f'TaDReP Database Creation v{tadrep.__version__}')
        print('Options and arguments:')
        print(f'\taction: {"create" if not update_path else "update"}')
        print(f'\tdatabase: {db_name}')
        print(f"\toutput: {output_path}")
        print(f'\ttmp directory: {tmp_path}')

    if(update_path or db_name == du.CUSTOM):
        if(not input_files):
            log.error('no files for database provided!')
            sys.exit('ERROR: No files for database provided!')
        input_files = [tu.check_file_permission(file, 'fasta') for file in input_files]

    fasta_tmp_path = tmp_path.joinpath(f'{db_name}.fna')
    log.info('merged fasta file: path=%s', fasta_tmp_path)

    if(update_path):
        db_tmp_path = update_path.joinpath('db')
        log.info('reverse database: path=%s', db_tmp_path)
        du.reverse_database(fasta_tmp_path, db_tmp_path, tmp_path)

    db_output_path = output_path.joinpath(f'{db_name}')
    if(db_output_path.exists()):
        if(force or update_path):
            shutil.rmtree(db_output_path)
            log.info('directory removed: path=%s', db_output_path)
        else:
            shutil.rmtree(str(tmp_path))
            log.debug('removed tmp dir: %s', tmp_path)
            log.debug('database directory already exists! path=%s', db_output_path)
            sys.exit(f'ERROR: Directory "{db_output_path.stem}" already exists in {output_path} !')

    db_output_path.mkdir(parents=True, exist_ok=True)
    log.info('directory created: path=%s', db_output_path)

    print('Database creation starting...')
    if(update_path):
        print('updating database...')
        with fasta_tmp_path.open('a') as fh_out:
            log.info('update file: path=%s', fasta_tmp_path)
            for file in input_files:
                log.info('add file: file=%s', file)
                with file.open('r') as fh_in:
                    fh_out.write(fh_in.read())
    else:
        if(selected_db == du.REFSEQ):
            dr.download_database(fasta_tmp_path)
        elif(selected_db == du.PLSDB):
            dp.download_database(fasta_tmp_path)
        else:
            print('Combining files...')
            with fasta_tmp_path.open('w+') as fh_out:
                for file in input_files:
                    log.info('write file: file=%s', file)
                    with file.open('r') as fh_in:
                        fh_out.write(fh_in.read())
    print('Create blast database...')
    db_file_name = db_output_path.joinpath('db')
    log.info('database files: name=%s', db_file_name)
    du.create_blast_db(db_file_name, fasta_tmp_path, tmp_path)
    
    
    tsv_output_path = db_output_path.joinpath('db.tsv')
    log.info('TSV file: name=%s, path=%s', tsv_output_path.stem, tsv_output_path)
    du.create_tsv(fasta_tmp_path, tsv_output_path)

    json_path = db_output_path.joinpath(f'{db_name}.json')
    log.info('JSON file: name=%s, path=%s', json_path.stem, json_path)
    db_plasmids = tio.import_sequences(fasta_tmp_path, sequence=True)
    plasmid_count = 1
    json_plasmids = {}

    for plasmid in db_plasmids.values():
        plasmid['file'] = db_name
        plasmid['old_id'] = plasmid['id']
        plasmid['id'] = str(plasmid_count)
        json_plasmids[plasmid_count] = plasmid

        plasmid_count += 1
    db_data = {'plasmids': json_plasmids}
    tio.export_json(db_data, json_path)
    
    print(f'Database successfully created\nDatabase path: {db_output_path}')

    shutil.rmtree(str(tmp_path))
    log.debug('removed tmp dir: %s', tmp_path)


if(__name__ == '__main__'):
    main()
