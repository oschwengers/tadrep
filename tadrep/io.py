import logging

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
    log.info('write genome sequences: path=%s, description=%s, wrap=%s', fasta_path, description, wrap)

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


def import_tsv(tsv_path):
    complete_path = Path(f'{tsv_path}/db.tsv')
    plasmids = {}
    with complete_path.open('r') as fh:
        for line in fh:
            cols = line.strip().split('\t')
            plasmid = {
                'id': cols[0],
                'description': cols[1],
                'length': int(cols[2])
            }
            plasmids[plasmid['id']] = plasmid
            log.info('imported: id=%s, description=%s, length=%s', plasmid['id'], plasmid['description'], plasmid['length'])
    return plasmids
