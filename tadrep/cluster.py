import logging

import tadrep.io as tio
import tadrep.config as cfg
import tadrep.utils as tu

log = logging.getLogger('CLUSTER')

def cluster_plasmids():
    # load json
    db_path = cfg.output_path.joinpath('db.json')
    db_data = tio.load_data(db_path)

    cluster = []
    # define method / criteria for clustering

    if(cfg.skip_cluster):
        # build clusters
        for plasmid in db_data['plasmids'].values():
            new_cluster = {
                'representative': plasmid['id'],
                'members': [plasmid['id']]
            }
            cluster.append(new_cluster)

        db_data['cluster'] = cluster
        prepare_blastdb(db_data)

    # find cluster representative

    # write json
    tio.export_json(db_data, db_path)


def prepare_blastdb(db_data):

    representative_ids = []

    # Get representatives from JSON file
    for cluster in db_data['cluster']:
        rep_id = cluster['representative']
        representative_ids.append(rep_id)

    reference_plasmids = {id: db_data['plasmids'][id] for id in representative_ids}

    cfg.verbose_print(f"Found {len(representative_ids)} representative plasmid(s)")
    log.info("Found %d representative plasmid(s)", len(representative_ids))

    # write multifasta for blast search
    fasta_tmp_path = cfg.tmp_path.joinpath("db.fasta")
    tio.export_sequences(reference_plasmids.values(), fasta_tmp_path)

    db_output_path = cfg.output_path.joinpath('reference_db/db')
    create_blast_db(db_output_path, fasta_tmp_path, cfg.tmp_path)


def create_blast_db(output_path, fasta_tmp_path, tmp_path):
    cmd_blastdb = [
        'makeblastdb',
        '-in', str(fasta_tmp_path),
        '-dbtype', 'nucl',
        '-out', str(output_path)
    ]
    log.debug('cmd=%s', cmd_blastdb)
    tu.run_cmd(cmd_blastdb, tmp_path)
