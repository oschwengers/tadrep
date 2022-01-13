import logging

import dnaplotlib as dpl
import matplotlib.pyplot as plt

from matplotlib import gridspec

logging.getLogger('matplotlib.font_manager').disabled = True
log = logging.getLogger('VISUALS')


HEIGHT = 0.5                                                                # Height per track
LABEL_Y_OFFSET = -3                                                         # Offset Label from track center
TRACK_EDGE_DISTANCE = 1                                                     # Distance between label and top of previous track
TRACK_CENTER_DISTANCE = 1.5 + abs(LABEL_Y_OFFSET) + TRACK_EDGE_DISTANCE     # Vertical distance between two plotted tracks
MARGIN_SPACE = 4                                                            # Vertical spacing between plotted elements and plot margin
LABEL_WRAPPING = 7                                                          # Number of characters per contig track
PLASMID_HEAD_MODIFIER = 100                                                 # Modifier of plasmid head length
CONTIG_HEAD_MODIFIER = 200                                                  # Modifier of contig head length


def create_plasmid_figure(plasmid, file_name, output_path):
    plasmid_head_length = int(plasmid['length'] / PLASMID_HEAD_MODIFIER)
    log.info('create graphic: plasmid=%s, head length=%i, output_path=%s', plasmid['id'], plasmid_head_length, output_path)

    plasmid_track = create_part_config(plasmid['id'], 0, plasmid['length'], label=True, head_length=plasmid_head_length)
    combined_tracks, num_contig_tracks = create_contig_tracks(sorted(plasmid['hits'], key=lambda k: k['reference_plasmid_start']), plasmid['length'])

    figure, grid = plot_setup(num_contig_tracks + 1, 2)
    plasmid_plot = plot_plasmid(grid[-1], [plasmid_track])
    label_width = plot_contigs(grid[0], plasmid_plot, num_contig_tracks, combined_tracks, y_label=file_name)
    finalize_plot(figure, output_path, label_width)


def create_contig_tracks(contig_hits, plasmid_length):

    head_length = int(plasmid_length / CONTIG_HEAD_MODIFIER)
    log.info('contig tracks: plasmid-length=%i, head-length=%i', plasmid_length, head_length)

    contig_complete_tracks = [[]]
    contig_matched_tracks = [[]]

    for hit in contig_hits:
        contig_fwd = True if hit['strand'] == '+' else False
        if(contig_fwd):
            start = hit['reference_plasmid_start'] - hit['contig_start'] + 1
            end = start + hit['contig_length']
            match_start = hit['reference_plasmid_start']
            match_end = hit['reference_plasmid_end']
        else:
            start = hit['reference_plasmid_end'] + hit['contig_start'] - 1
            end = start - hit['contig_length']
            match_start = hit['reference_plasmid_end']
            match_end = hit['reference_plasmid_start']

        if(abs(end - start) < head_length):
            log.info('contig skipped: id=%s, length=%s', hit['contig_id'], hit['contig_length'])
            continue

        contig_config = create_part_config(hit['contig_id'], int(start), int(end), fwd=contig_fwd, head_length=head_length)
        match_config = create_part_config(hit['contig_id'], int(match_start), int(match_end), fwd=contig_fwd, color=(0.38, 0.82, 0.32), label=True, label_y_offset=LABEL_Y_OFFSET, head_length=head_length)

        track = find_track(contig_complete_tracks, contig_matched_tracks, contig_config, len(contig_hits) + 1)
        contig_complete_tracks[track].append(contig_config)
        contig_matched_tracks[track].append(match_config)
        log.info('contig added: id=%s, start=%s, end=%s, length=%s, type=%s, strand=%s, fwd=%s',
                  contig_config['name'], contig_config['start'], contig_config['end'], hit['length'], contig_config['type'], hit['strand'], contig_config['fwd'])
    combined_tracks = combine_tracks(contig_complete_tracks, contig_matched_tracks)

    return combined_tracks, len(contig_complete_tracks)


def plot_setup(num_total_tracks, num_grids):
    height = HEIGHT * num_total_tracks
    fig = plt.figure(figsize=(10.0, height))
    gs = gridspec.GridSpec(num_grids, 1, height_ratios=[(num_total_tracks - 1) * HEIGHT, HEIGHT + 0.1], figure=fig)
    gs.update(hspace=0.2)
    log.info('figure created: tracks=%i, height=%i, grids=%i', num_total_tracks, height, num_grids)
    return fig, gs


def plot_plasmid(grid_pos, design):
    plasmid_plot = plt.subplot(grid_pos)
    dr = dpl.DNARenderer(linewidth=0.8)
    part_renderers = dr.trace_part_renderers()
    start, end = dr.renderDNA(plasmid_plot, design, part_renderers, plot_backbone=False)
    plasmid_plot.set_xlim([start - 100, end + 100])
    plasmid_plot.set_ylim([-(MARGIN_SPACE - 1), MARGIN_SPACE - 1])
    plasmid_plot.axis('off')
    log.info('subplot created: plasmid, start=%s, end=%s, y_limitation=%s', start, end, MARGIN_SPACE - 1)
    return plasmid_plot


def plot_contigs(grid_pos, plasmid_plot, num_contig_tracks, combined_tracks, y_label=None):
    contig_plot = plt.subplot(grid_pos, sharex=plasmid_plot)
    set_ylim_contigs = TRACK_CENTER_DISTANCE * (num_contig_tracks - 1) + MARGIN_SPACE
    dr = dpl.DNARenderer(linewidth=0.8)
    part_renderers = dr.trace_part_renderers()

    start, end = dr.renderDNA(contig_plot, combined_tracks, part_renderers, plot_backbone=False)
    contig_plot.set_ylim([-MARGIN_SPACE, set_ylim_contigs])  # set view limitations to show everything
    log.info('subplot created: contig, start=%s, end=%s, tracks=%s', start, end, num_contig_tracks)
    if(y_label):  # remove ticks and axis
        wrapped_label, label_width = wrap_y_label(y_label, num_contig_tracks)
        contig_plot.set_ylabel(wrapped_label, labelpad=1, fontsize=8)
        contig_plot.spines['right'].set_visible(False)
        contig_plot.spines['top'].set_visible(False)
        contig_plot.spines['bottom'].set_visible(False)
        contig_plot.spines['left'].set_visible(False)
        contig_plot.set_xticks([])
        contig_plot.set_yticks([])
        log.debug('Y_label set: pos=%s, label=%s', grid_pos, y_label)
    else:
        label_width = 0
        contig_plot.axis('off')
        log.debug('subplot axis: disabled')
    return label_width


def finalize_plot(figure, output_path, pad_left=0):
    plt.subplots_adjust(left=pad_left + 0.01, right=0.99, top=0.985, bottom=0.01)
    figure.savefig(output_path, dpi=600, format='pdf')
    plt.close('all')
    log.info('plot exported: path=%s', output_path)


def wrap_y_label(y_label, num_contig_tracks):
    wrapping = LABEL_WRAPPING * num_contig_tracks
    lines = []
    for i in range(0, len(y_label), wrapping):
        lines.append(y_label[i:i + wrapping])
    label_width = len(lines) * 0.01
    return '\n'.join(lines), label_width


def create_part_config(name, start, end, head_length=10, label_y_offset=0, color=(0.70, 0.70, 0.70), fwd=True, label=False):
    part_config = {
        'type': 'CDS',
        'name': name,
        'start': start,
        'end': end,
        'fwd': fwd,
        'opts': {
            'label': name if label else None,
            'arrowhead_height': 0,
            'arrowhead_length': head_length,
            'label_y_offset': label_y_offset,
            'label_style': 'italic',
            'color': color
        }
    }
    log.debug('design created! id=%s, start=%s, end=%s', part_config['name'], part_config['start'], part_config['end'])
    return part_config


def combine_tracks(contig_hits, match_hits):
    combined_tracks = []
    for count, track in enumerate(contig_hits):
        for part in track:
            part['opts']['y_offset'] = count * TRACK_CENTER_DISTANCE
        combined_tracks = combined_tracks + track
    for count, track in enumerate(match_hits):
        for part in track:
            part['opts']['y_offset'] = count * TRACK_CENTER_DISTANCE
        combined_tracks = combined_tracks + track
    return combined_tracks


def find_track(contig_track, match_track, contig_config, length):
    for track in range(length + 1):
        if(not contig_track[track]):
            return track

        previous_contig = contig_track[track][-1]
        if(previous_contig['fwd']):
            if(previous_contig['end'] < contig_config['start'] and previous_contig['end'] < contig_config['end']):
                return track
        else:
            if(previous_contig['start'] < contig_config['start'] and previous_contig['start'] < contig_config['end']):
                return track

        if(track + 1 == len(contig_track)):
            contig_track.append([])
            match_track.append([])
