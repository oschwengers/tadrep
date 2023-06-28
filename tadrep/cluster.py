import logging

import tadrep.io as tio
import tadrep.config as cfg
import tadrep.utils as tu

log = logging.getLogger('CLUSTER')

def cluster_plasmids():
    # load json
    db_path = cfg.output_path.joinpath('db.json')
    db_data = tio.load_data(db_path)


    # write multifasta
    fasta_path = cfg.tmp_path.joinpath('plasmids.extracted.fna')
    tio.export_sequences(db_data['plasmids'].values(), fasta_path)

    # cluster sequences
    cmd_cdhitest = [
        'cd-hit-est',
        '-i', str(fasta_path),
        '-o', str(cfg.tmp_path.joinpath('plasmids.clustered')),
        '-c', '0.8',
        '-S', '1000',
        '-g', '1'
    ]
    log.debug('cmd=%s', cmd_cdhitest)
    tu.run_cmd(cmd_cdhitest, cfg.tmp_path)

    # parse clusters
    clusters = []
    current_cluster = None
    with cfg.tmp_path.joinpath('plasmids.clustered.clstr').open('r') as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if current_cluster is not None:
                    clusters.append(current_cluster)
                current_cluster = {
                    'id': f"p{line[1:].split(' ')[1]}",  # cluster number
                    'representative': None,
                    'members': []
                }
            else:
                plasmid_id = line.split('\t')[1].split(' ')[1][1:-3]  # exmaple: "0	6222nt, >And12566.short-33... *"
                current_cluster['members'].append(plasmid_id)
                if '*' in line:
                    current_cluster['representative'] = plasmid_id
        if current_cluster is not None:
            clusters.append(current_cluster)

    db_data['clusters'] = clusters

    cfg.verbose_print('Detected plasmid clusters:')
    for cluster in clusters:
        members = ', '.join(cluster['members'])
        cfg.verbose_print(f"Cluster {cluster['id']}\n\tsize: {len(cluster['members'])}\n\trepresentative: {cluster['representative']}\n\tmembers: {members}")
        log.info('cluster: %s, size: %d, rep: %s, members: %s', cluster['id'], len(cluster['members']), cluster['representative'], cluster['members'])
    
    # write json
    tio.export_json(db_data, db_path)