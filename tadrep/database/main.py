import logging
import shutil
import sys

import tadrep.io as tio
import tadrep.config as cfg
import tadrep.database.refseq as dr
import tadrep.database.plsdb as dp


log = logging.getLogger('DATABASE')


def create_database():
    cfg.verbose_print(f'\tdatabase: {cfg.db_type}')
    cfg.verbose_print(f"\toutput: {cfg.output_path}")
    cfg.verbose_print(f'\ttmp directory: {cfg.tmp_path}')

    fasta_tmp_path = cfg.tmp_path.joinpath(f'{cfg.db_type}.fna')
    log.info('merged fasta file: path=%s', fasta_tmp_path)

    db_output_path = cfg.output_path.joinpath(f'{cfg.db_type}')
    if(db_output_path.exists()):
        if(cfg.force):
            shutil.rmtree(db_output_path)
            log.info('directory removed: path=%s', db_output_path)
        else:
            shutil.rmtree(str(cfg.tmp_path))
            log.debug('removed tmp dir: %s', cfg.tmp_path)
            log.debug('database directory already exists! path=%s', db_output_path)
            sys.exit(f'ERROR: Directory "{db_output_path.stem}" already exists in {cfg.output_path} !')

    db_output_path.mkdir(parents=True, exist_ok=True)
    log.info('directory created: path=%s', db_output_path)

    print('Database creation starting...')
    if(cfg.db_type == 'refseq'):
        dr.download_database(fasta_tmp_path)
    else:
        dp.download_database(fasta_tmp_path)

    print('Create JSON database...')
    json_path = db_output_path.joinpath(f'{cfg.db_type}.json')
    log.info('JSON database: name=%s, path=%s', json_path.stem, json_path)
    db_plasmids = tio.import_sequences(fasta_tmp_path, sequence=True)
    json_plasmids = {}

    for plasmid in db_plasmids.values():
        plasmid['file'] = cfg.db_type
        json_plasmids[plasmid['id']] = plasmid

    db_data = {
        'plasmids': json_plasmids
        }
    tio.export_json(db_data, json_path)

    print(f'Database successfully created\nDatabase path: {db_output_path}')
