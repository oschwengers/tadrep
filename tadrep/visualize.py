import logging

from pygenomeviz import GenomeViz

import tadrep.config as cfg
import tadrep.io as tio

logging.getLogger('matplotlib.font_manager').disabled = True
log = logging.getLogger('VISUALIZE')


def plot():
    db_path = cfg.output_path.joinpath('db.json')
    db_data = tio.load_data(db_path)

    for cluster in db_data['clusters']:
        detected_genomes = cluster.get('found_in', {})
        if(len(detected_genomes) == 0):
            continue
        for draft_genome, hits in cluster['found_in'].items():
            cfg.verbose_print(f"Plasmid: {cluster['id']}, genome: {draft_genome}, hits: {len(hits)}")
            log.info('plasmid: %s, genome: %s, contig-hits: %d', cluster['id'], draft_genome, len(hits))
            output_path = cfg.output_path.joinpath(f"{draft_genome}-{cluster['id']}.pdf")
            create_figure(cluster['id'], db_data['plasmids'][cluster['representative']]['length'], hits, output_path)


def create_figure(plasmid_id, plasmid_length, hits, output_path):
    gv = GenomeViz(tick_style='axis')
    track = gv.add_feature_track(plasmid_id, plasmid_length)
    min_size = int(plasmid_length * (cfg.omit_ratio / 100))
    log.info('Skipping contigs shorter than: %s', min_size)

    for hit in hits:
        if(hit['length'] < min_size):
            log.info('Skipped contig %s, length: %d', hit['contig_id'], hit['length'])
            continue
        strand = 1 if hit['strand'] == '+' else -1
        contig_colour = get_gradient_colour(hit['perc_identity'])
        log.debug('contig: %s, plasmid_start: %d, plasmid_end: %d, identity: %f, colour: %f', hit['contig_id'], hit['reference_plasmid_start'], hit['reference_plasmid_end'], hit['perc_identity'], contig_colour)
        track.add_feature(hit['reference_plasmid_start'], hit['reference_plasmid_end'], strand, label=hit['contig_id'],
         plotstyle=cfg.plot_style, labelsize=cfg.label_size, labelcolor=cfg.label_color, facecolor=str(contig_colour), linewidth=cfg.line_width, 
         labelrotation=cfg.label_rotation, labelhpos=cfg.label_hpos, labelha=cfg.label_ha, arrow_shaft_ratio=cfg.arrow_shaft_ratio, size_ratio=cfg.size_ratio)

    fig = gv.plotfig()
    fig.savefig(output_path, dpi=600, format='pdf')


def get_gradient_colour(plasmid_coverage):
    colour = 1.0
    for interval in range(cfg.interval_number + 1):
        if(1.0 - interval * cfg.interval_size <= plasmid_coverage):
            colour = interval / cfg.interval_number
            break
    return colour
