import logging

from pygenomeviz import GenomeViz

import tadrep.config as cfg
import tadrep.io as tio

log = logging.getLogger('VISUALIZE')


def plot():

    db_path = cfg.output_path.joinpath('plasmids.json')
    db_data = tio.load_data(db_path)

    for plasmid in db_data.values():
        for draft_genome in plasmid["found_in"].keys():
            hits = plasmid['found_in'][draft_genome]
            cfg.verboseprint(f'Plasmid ID: {plasmid["id"]}, Draft genome: {draft_genome}, Hits: {len(hits)}')
            log.info('plasmid: %s, genome: %s, contig-hits: %d', plasmid['id'], draft_genome, len(hits))
            output_path = cfg.output_path.joinpath(f'{draft_genome}-{plasmid["id"]}.pdf')
            create_figure(plasmid['id'], plasmid['length'], hits, output_path)


def create_figure(plasmid_id, plasmid_length, hits, output_path):
    gv = GenomeViz(tick_style='axis')
    track = gv.add_feature_track(plasmid_id, plasmid_length)
    for hit in hits:
        strand = 1 if hit['strand'] == '+' else -1
        track.add_feature(hit['reference_plasmid_start'], hit['reference_plasmid_end'], strand, label=hit['contig_id'], plotstyle='arrow')

    fig = gv.plotfig()
    fig.savefig(output_path, dpi=600, format='pdf')
