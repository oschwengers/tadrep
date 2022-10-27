import logging
import json

from pathlib import Path

from Bio import SeqIO
from xopen import xopen


log = logging.getLogger('IO')


FASTA_LINE_WRAPPING = 60


def import_sequences(contigs_path, sequence=False):
    """Import sequences."""
    contigs = {}
    with xopen(str(contigs_path), threads=0) as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            seq = str(record.seq).upper()
            contig = {
                'id': record.id,
                'description': record.description.split(' ', maxsplit=1)[1] if ' ' in record.description else '',
                'sequence': seq if sequence else None,
                'length': len(seq)
            }
            log.info(
                'imported: id=%s, length=%i, description=%s',
                contig['id'], contig['length'], contig['description']
            )
            contigs[contig['id']] = contig
    return contigs


def export_sequences(contigs, fasta_path, description=False, wrap=False):
    """Write sequences to Fasta file."""
    log.debug('write: path=%s, description=%s, wrap=%s', fasta_path, description, wrap)
    with fasta_path.open('wt') as fh:
        for contig in contigs:
            if(description):
                fh.write(f">{contig['id']} {contig['description']}\n")
            else:
                fh.write(f">{contig['id']}\n")
            if(wrap):
                fh.write(wrap_sequence(contig['sequence']))
            else:
                fh.write(contig['sequence'])
                fh.write('\n')


def wrap_sequence(sequence):
    lines = []
    for i in range(0, len(sequence), FASTA_LINE_WRAPPING):
        lines.append(sequence[i:i + FASTA_LINE_WRAPPING])
    return '\n'.join(lines) + '\n'


def import_tsv(database_path):
    plasmids = {}
    complete_path = database_path.joinpath('db.tsv')
    with complete_path.open('r') as fh:
        for line in fh:
            cols = line.strip().split('\t')
            plasmid = {
                'id': cols[0],
                'description': cols[1],
                'length': int(cols[2])
            }
            plasmids[plasmid['id']] = plasmid
            log.debug('imported: id=%s, description=%s, length=%i', plasmid['id'], plasmid['description'], plasmid['length'])
    log.info('imported: # plasmids=%i', len(plasmids))
    return plasmids


def import_json(json_path):
    with open(json_path, 'r') as fh:
        data = json.load(fh)
    log.info('imported json: # sequences=%i', len(data))
    return data


def export_json(data, json_path):
    log.info('write json: path=%s, # sequences=%i', json_path, len(data))
    with open(json_path, 'w') as fh:
        json.dump(data, fh)


def check_existing_JSON(json_path):
    if json_path.is_file():
        log.info('%s existing', json_path)
        plasmid_dict = import_json(json_path)
    else:
        log.info('%s NOT existing', json_path)
        plasmid_dict = {}
    return plasmid_dict
