import tadrep.config as cfg
import tadrep.io as tio
import tadrep.visuals as tv
import tadrep.blast as tb
import tadrep.plasmids as tp

import numpy as np
import logging
import sys
import multiprocessing as mp

log = logging.getLogger('DETECTION')
verboseprint = print if cfg.verbose else lambda *a, **k: None

def detect_and_reconstruct():

    ############################################################################
    # Import plasmid sequences
    # - parse contigs in Fasta file
    ############################################################################
    if(cfg.plasmids_path):
        try:
            verboseprint('\nimport plasmids sequences...')
            reference_plasmids = tio.import_sequences(cfg.plasmids_path, sequence=True)
            log.info('imported reference plasmids: sequences=%i, file=%s', len(reference_plasmids), cfg.plasmids_path)
            for ref_plasmid in sorted(reference_plasmids.values(), key=lambda k: k['length']):
                size = ref_plasmid['length']
                seq = ref_plasmid['sequence']
                gc = (seq.count('G') + seq.count('C')) / (size - seq.count('N'))
                verboseprint(f"\t{ref_plasmid['id']}: {size:,} bp, {(gc*100):0.1f} % GC")
            verboseprint(f'\n\ttotal: {len(reference_plasmids)}')
        except ValueError:
            log.error('wrong reference plasmids file format!', exc_info=True)
            sys.exit('ERROR: wrong reference plasmids file format!')
    else:
        try:
            verboseprint('\nload reference plasmids database...')
            reference_plasmids = tio.import_tsv(cfg.database_path)
            log.info('imported reference plasmids: sequence=%i, file=%s', len(reference_plasmids), cfg.database_path)
            verboseprint(f'\timported: {len(reference_plasmids)}')
        except:
            log.error('wrong database path!', exc_info=True)
            sys.exit('ERROR: wrong database path!')
    

    ############################################################################
    # Prepare summary output file
    # - create file
    # - write header into file
    ############################################################################
    plasmid_dict = {}
    plasmid_string_summary = []

    verboseprint('Analyze genome sequences...')
    values = ((genome_path, reference_plasmids, genome_index) for genome_index, genome_path in enumerate(cfg.genome_path))
    with mp.Pool(cfg.threads) as pool:
        genomes_summary = pool.starmap(pooling, values)

    for genome_index, plasmid_summary in genomes_summary:
        for plasmid_id, summary_string in plasmid_summary.items():
            plasmid_string_summary.append(summary_string)
            if(plasmid_id not in plasmid_dict):
                plasmid_dict[plasmid_id] = [0] * len(cfg.genome_path)
                log.info('Plasmid added: id=%s', plasmid_id)
            plasmid_dict[plasmid_id][genome_index] = 1

    cfg.summary_path = cfg.output_path.joinpath('summary.tsv')
    with cfg.summary_path.open('w') as fh:
        fh.write(f'# {len(cfg.genome_path)} draft genome(s), {len(reference_plasmids)} reference plasmid(s)\n')
        fh.write('Genome\tPlasmid\tCoverage\tIdentity\tContigs\tContig IDs\n')
        for line in plasmid_string_summary:
            fh.write(line)

    if(plasmid_dict):
        write_cohort_table(plasmid_dict)


def pooling(genome, reference_plasmids, index):
    log = logging.getLogger('PROCESS')

    # Import draft genome contigs
    try:
        contigs = tio.import_sequences(genome, sequence=True)
        log.info('imported genome contigs: genome=%s, # contigs=%i', genome, len(contigs))
    except ValueError:
        log.error('wrong genome file format!', exc_info=True)
        sys.exit('ERROR: wrong genome file format!')

    sample = genome.stem
    blast_output_path = cfg.tmp_path.joinpath(f'{sample}-{index}-blastn.tsv')
    hits = tb.search_contigs(genome, blast_output_path)  # plasmid raw hits
    filtered_hits = tb.filter_contig_hits(sample, hits, reference_plasmids)  # plasmid hits filtered by coverage and identity
    detected_plasmids = tp.detect_reference_plasmids(sample, filtered_hits, reference_plasmids)  # detect reference plasmids above cov/id thresholds

    # Write output files
    plasmid_summary_strings = {}
    sample_summary_path = cfg.output_path.joinpath(f'{sample}-summary.tsv')
    with sample_summary_path.open('w') as ssp:
        ssp.write(f"plasmid\tcontig\tcontig start\tcontig end\tcontig length\tcoverage[%]\tidentity[%]\talignment length\tstrand\tplasmid start\tplasmid end\tplasmid length\n")

        for plasmid in detected_plasmids:

            plasmid_contigs_sorted = tp.reconstruct_plasmid(plasmid, contigs)
            
            prefix = f"{cfg.prefix}-{sample}-{plasmid['reference']}" if cfg.prefix else f"{sample}-{plasmid['reference']}"
            plasmid_contigs_path = cfg.output_path.joinpath(f'{prefix}-contigs.fna')
            plasmid_pseudosequence_path = cfg.output_path.joinpath(f'{prefix}-pseudo.fna')
            tio.export_sequences(plasmid_contigs_sorted, plasmid_contigs_path, description=True, wrap=True)
            tio.export_sequences([plasmid], plasmid_pseudosequence_path, description=True, wrap=True)

            # Create string for plasmid summary
            plasmid_summary_strings[plasmid['reference']] = f"{sample}\t{plasmid['reference']}\t{plasmid['coverage']:.3f}\t{plasmid['identity']:.3f}\t{len(plasmid['hits'])}\t{','.join([hit['contig_id'] for hit in plasmid['hits']])}\n"

            # Write detailed plasmid hits to sample summary file
            for hit in plasmid['hits']:
                ssp.write(f"{plasmid['reference']}\t{hit['contig_id']}\t{hit['contig_start']}\t{hit['contig_end']}\t{hit['contig_length']}\t{hit['coverage']:.3f}\t{hit['perc_identity']:.3f}\t{hit['length']}\t{hit['strand']}\t{hit['reference_plasmid_start']}\t{hit['reference_plasmid_end']}\t{plasmid['length']}\n")

            png_path = cfg.output_path.joinpath(f"{sample}-{plasmid['reference']}.pdf")
            tv.create_plasmid_figure(plasmid, genome.name, png_path)

    cfg.lock.acquire()
    log.debug('lock acquired: genome=%s, index=%s', sample, index)
    print(f'\n\nGenome: {sample}, contigs: {len(contigs)}, detected plasmids: {len(detected_plasmids)}')
    for plasmid in detected_plasmids:
        # Create console output
        if (cfg.verbose):
            print(f"\n\tplasmid: {plasmid['reference']}, length: {plasmid['length']} bp, contig hits: {len(plasmid['hits'])}, coverage: {plasmid['coverage'] * 100:1.1f}%, identity: {plasmid['identity'] * 100:1.1f}%")
            print(f"\t{'contig':^17}  alignment length  contig length  contig start  contig end  strand  plasmid start  plasmid end  coverage[%]  identity[%]")
            for hit in plasmid['hits']:
                contig = contigs[hit['contig_id']]
                print(f"\t{hit['contig_id']:^18} {hit['length']:>11} {contig['length']:>14} {hit['contig_start']:>13} {hit['contig_end']:>13} {hit['strand']:^13} {hit['reference_plasmid_start']:>9} {hit['reference_plasmid_end']:>12} {(hit['length'] / contig['length'] * 100):>13.1f} {hit['perc_identity'] * 100:>12.1f}")
        else:
            print(f"{sample}\t{plasmid['id']}\t{plasmid['length']}\t{len(plasmid['hits'])}\t{plasmid['coverage']:f}\t{plasmid['identity']:f}\t{','.join(contig['contig_id'] for contig in plasmid['hits'])}")
    log.debug('lock released: genome=%s, index=%s', sample, index)
    cfg.lock.release()

    return index, plasmid_summary_strings


def write_cohort_table(plasmid_dict):
    plasmid_cohort_path = cfg.output_path.joinpath('plasmids.tsv')
    with plasmid_cohort_path.open('w') as fh:
        plasmid_order = []  # write plasmid header
        for plasmid in plasmid_dict.keys():
            plasmid_order.append(plasmid_dict[plasmid])
            fh.write(f'\t{plasmid}')
        fh.write('\n')
        
        transposed_plasmid_order = np.array(plasmid_order).T.tolist()  # transpose information for easier writing
        for num_genome, genome in enumerate(cfg.genome_path):  # mark which plasmid was found for each draft genome
            sample = genome.stem
            fh.write(f'{sample}')
            for plasmid in transposed_plasmid_order[num_genome]:
                fh.write(f'\t{"X" if plasmid else "-"}')
            fh.write('\n')
