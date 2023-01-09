import logging

import tadrep.io as tio
import tadrep.config as cfg

log = logging.getLogger('CLUSTER')

def cluster():
    # load json
    db_path = cfg.output_path.joinpath('db.json')
    db_data = tio.load_data(db_path)

    plasmids = db_data.get('plasmids', {})
    cluster = db_data.get('cluster', [])
    # define method / criteria for clustering

    # build clusters
    for plasmid in plasmids.values():
        new_cluster = {
            'representative': plasmid['new_id'],
            'members': [plasmid['new_id']]
        }
        cluster.append(new_cluster)

    db_data['cluster'] = cluster
    # find cluster representative
    
    # write json
    tio.export_json(db_data, db_path)
