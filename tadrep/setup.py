import logging
import re

import tadrep.config as cfg
import tadrep.utils as tu

log = logging.getLogger('SETUP')

description_pattern = re.compile('^>([a-zA-Z0-9()]+)_([0-9]+)_(.*)_(.*)$')

def download_inc_reference():

    inc_types_path = cfg.output_path.joinpath('inc-types.fasta')

    # check if inc-types.fasta is already available
    if(inc_types_path.is_file()):
        print(f'{inc_types_path.name} already exists in {cfg.output_path}')
        log.debug('Found inc-types.fasta in path %s', cfg.output_path)
        return

    plasmidfinder_git = 'https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git'

    print(f'Downloading {inc_types_path.name}')
    log.debug('Missing %s, cloning from %s', inc_types_path.name, plasmidfinder_git)

    git_cmd = ['git', 'clone', str(plasmidfinder_git)]
    tu.run_cmd(git_cmd, cfg.tmp_path)

    with open(inc_types_path, 'w') as inc_types:

        # read each .fsa file
        plasmidfinder_path = cfg.tmp_path.joinpath('plasmidfinder_db')
        for fasta_file in plasmidfinder_path.glob('*.fsa'):
            log.debug('Reading %s', fasta_file)
            with open(fasta_file, 'r') as inc_raw:

                for line in inc_raw.readlines():
                    # if line is header
                    if(line.startswith('>')):
                        # check if pattern
                        if(description_pattern.match(line)):
                            inc_types.write(line.split('_')[0])
                        else:
                            inc_types.write(line)
                    # write seq in uppercase
                    else:
                        inc_types.write(line.upper())
                    inc_types.write('\n')
