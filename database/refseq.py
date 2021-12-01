import logging
import sys
import urllib.error
import urllib.request
import gzip


log = logging.getLogger('REFSEQ')


NCBI_PATH = 'https://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid'
FILE_NUMBERS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]


def download_database(merged_file_path):
    log.info('Download NCBI files: destination=%s', merged_file_path)
    with merged_file_path.open('wb') as fh_out:
        for file in FILE_NUMBERS:
            ncbi_url = f'{NCBI_PATH}/plasmid.{file}.1.genomic.fna.gz'
            print(f'Downloading file: {ncbi_url} ...')
            log.info('Download file: number=%s, url=%s, destination=%s', file, ncbi_url, merged_file_path)
            try:
                with urllib.request.urlopen(ncbi_url) as response:
                    with gzip.GzipFile(fileobj=response) as fh_in:
                        fh_out.write(fh_in.read())
            except urllib.error.URLError:
                log.debug('URLError occurred! Could not read %s', ncbi_url)
                sys.exit(f'ERROR: Could not read {ncbi_url}')
            except OSError:
                log.error('Could not write file %s', merged_file_path)
                sys.exit(f'ERROR: Could not write {merged_file_path}')
