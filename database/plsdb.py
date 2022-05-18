import logging
import sys
import urllib.request
import urllib.error
import bz2
import shutil


log = logging.getLogger('PLSDB')


PLSDB_URL = 'https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/plsdb.fna.bz2'
FILE = 'plsdb.fna'


def download_database(tmp_path):
    log.info('download PLSDB files')
    print(f'Downloading file: {PLSDB_URL} ...')
    try:
        with urllib.request.urlopen(PLSDB_URL) as url_in:
            with bz2.BZ2File(url_in) as fh_in:
                with open(tmp_path, "wb") as fh_out:
                    log.info('decompressing file: destination=%s', tmp_path)
                    shutil.copyfileobj(fh_in, fh_out)
        log.info('copying file finished')
    except urllib.error.URLError:
        log.debug('URLError occurred! Could not resolve %s, error=%s', PLSDB_URL)
        sys.exit(f'ERROR: Could not resolve {PLSDB_URL}')
    except OSError:
        log.error('could not write file')
        sys.exit(f'ERROR: Could not write file')
