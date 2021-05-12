import logging
import sys
import urllib.request
import urllib.error
import zipfile
import io

import tadrep.utils as tu

log = logging.getLogger('PLSDB')

PLSDB_URL = 'https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/?zip'
FILES = ['plsdb.fna.nhr', 'plsdb.fna.nin', 'plsdb.fna.nsq']


def download_database(tmp_path):
    log.info('Downloading PLSDB files')
    try:
        with urllib.request.urlopen(PLSDB_URL) as url_in:
            with zipfile.ZipFile(io.BytesIO(url_in.read())) as fh_in:
                for file in FILES:
                    output_file = tmp_path.joinpath(f'{file}')
                    with output_file.open('wb') as fh_out:
                        fh_out.write(fh_in.read(file))
                    log.info('Fetching file: file=%s, destination=%s', file, output_file)
        log.info('Copying files finished')
    except urllib.error.URLError:
        log.debug('URLError occurred! Could not resolve %s, error=%s', PLSDB_URL)
        sys.exit(f'ERROR: Could not resolve {PLSDB_URL}')
    except OSError:
        log.error('Could not write file')
        sys.exit(f'ERROR: Could not write file')


def reverse_database(fasta_path, tmp_path):
    database_path = tmp_path.joinpath('plsdb.fna')
    log.debug('Database path: input=%s, output=%s', database_path, fasta_path)

    cmd_reverse = [
        'blastdbcmd',
        '-entry', 'all',
        '-db', str(database_path),
        '-out', str(fasta_path)
    ]
    log.debug('cmd_reverse=%s', cmd_reverse)
    tu.run_cmd(cmd_reverse, tmp_path)
