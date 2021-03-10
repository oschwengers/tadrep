
import logging

from Bio import SeqIO
from xopen import xopen

import tadrep.constants as bc


log = logging.getLogger('FASTA')

FASTA_LINE_WRAPPING = 60


def import_sequences(contigs_path):
    """Import sequences."""
    contigs = []
    with xopen(str(contigs_path), threads=0) as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            seq = str(record.seq).upper()
            contig = {
                'id': record.id,
                'description': record.description,
                'sequence': seq,
                'length': len(seq)
            }
            log.info(
                'imported: id=%s, length=%i, description=%s',
                contig['id'], contig['length'], contig['description']
            )
            contigs.append(contig)
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
