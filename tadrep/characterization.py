import logging
import sys

import tadrep.io as tio
import tadrep.config as cfg
import tadrep.utils as tu


log = logging.getLogger('CHARACTERIZE')


def characterize():
    # load json from output path
    json_path = cfg.output_path.joinpath('db.json')
    existing_data = tio.load_data(json_path)
    plasmids = existing_data.get('plasmids', {})

    # check if data is available
    if(not existing_data or not plasmids):
        log.error('Failed to load data from %s', json_path)
        sys.exit(f'ERROR: Failed to load data from {json_path}!')

    for plasmid in plasmids:
        # set length
        plasmid['length'] = len(plasmid['sequence'])

        # set gc_content
        plasmid['gc_content'] = gc_content(plasmid['sequence'])

        # calculate INC_types (Platon)
            # write multifasta
            # find inc-types.fasta
            # threads = cfg.threads

        # gene prediction (pyrodigal)

    download_inc_types()
    # update json
    tio.export_json(existing_data, json_path)


def gc_content(sequence):
    # count C and G occurenxes
    c_content = sequence.count('C')
    g_content = sequence.count('G')
    gc_combined = c_content + g_content

    # calculate gc percentage
    percentage = gc_combined / float(len(sequence))

    return percentage


def download_inc_types():
    
    inc_types_path = cfg.output_path.joinpath('inc-types.fasta')
    
    # check if inc-types.fasta is already available
    if(inc_types_path.is_file()):
        return
    
    inc_types_cmd = [
        'git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git;',
        'for f in plasmidfinder_db/*.fsa; do cat $f >> inc-types-raw.fasta; done;',
        'cat inc-types-raw.fasta | sed -r "s/^>([a-zA-Z0-9()]+)_([0-9]+)_(.*)_(.*)$/>\1 \4/" > inc-types-mixed.fasta;',
        "awk '{if(!/>/){print toupper($0)}else{print $1}}' inc-types-mixed.fasta > ", f"{inc_types_path};",
        'rm -rf plasmidfinder_db inc-types-raw.fasta inc-types-mixed.fasta'
    ]

    tu.run_cmd(inc_types_cmd, cfg.output_path)
    return


def search_inc_types():
    pass


def gene_prediction():
    pass
