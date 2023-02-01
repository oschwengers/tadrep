import tadrep.config as cfg
import tadrep.io as tio
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
    # - write multi Fasta file
    ############################################################################

    if(not cfg.db_data.get('plasmids', {})):
        log.debug("No plasmids in %s !", cfg.db_path)
        sys.exit(f"ERROR: No plasmids in database {cfg.db_path}!")

    if(not cfg.db_data.get('cluster', [])):
        log.debug("No Clusters in %s!", cfg.db_path)
        sys.exit(f"ERROR: No cluster in database {cfg.db_path}")

    log.info("Loaded %d cluster with %d plasmids", len(cfg.db_data['cluster']), len(cfg.db_data['plasmids'].keys()))
    cfg.verboseprint("Loaded data:")
    cfg.verboseprint(f"\t{len(cfg.db_data['cluster'])} cluster")
    cfg.verboseprint(f"\t{len(cfg.db_data['plasmids'].keys())} plasmids total")

    representative_ids = []
    lookup_plasmid_ids = {}

    # Get representatives from JSON file
    for cluster in cfg.db_data['cluster']:
        rep_id = cluster['representative']
        representative_ids.append(rep_id)
        if(cfg.blastdb_path):
            lookup_plasmid_ids[cfg.db_data['plasmids'][rep_id]['old_id']] = rep_id      # lookuptable for old_id(plasmid id from file) and new id (given by tadrep)

    reference_plasmids = {id: cfg.db_data['plasmids'][id] for id in representative_ids}

    cfg.verboseprint(f"Found {len(representative_ids)} representative plasmid(s)")
    log.info("Found %d representative plasmid(s)", len(representative_ids))

    if(not cfg.blastdb_path):
        # write multifasta for blast search
        fasta_path = cfg.output_path.joinpath("db.fasta")
        tio.export_sequences(reference_plasmids.values(), fasta_path)

    ############################################################################
    # Prepare summary output file
    # - create file
    # - write header into file
    ############################################################################
    plasmid_dict = {}
    plasmid_string_summary = []
    plasmid_detected = {}

    cfg.verboseprint('Analyze genome sequences...')
    values = ((genome_path, reference_plasmids, genome_index, lookup_plasmid_ids) for genome_index, genome_path in enumerate(cfg.genome_path))
    with mp.Pool(cfg.threads) as pool:
        genomes_summary = pool.starmap(pooling, values)

    for genome_index, plasmid_summary in genomes_summary:
        for plasmid in plasmid_summary:
            plasmid_id = plasmid['reference']
            if(plasmid_id not in plasmid_detected):
                plasmid_detected[plasmid_id] = {k: plasmid[k] for k in ['id', 'reference', 'length']}
                plasmid_detected[plasmid_id]['found_in'] = {}
            plasmid_detected[plasmid_id]['found_in'][plasmid['genome']] = plasmid['hits']
            plasmid_detected[plasmid_id]['old_id'] = cfg.db_data['plasmids'][plasmid_id]['old_id']

            # Create string for plasmid summary
            plasmid_string_summary.append(f"{plasmid['genome']}\t{plasmid['reference']}\t{plasmid['coverage']:.3f}\t{plasmid['identity']:.3f}\t{len(plasmid['hits'])}\t{','.join([hit['contig_id'] for hit in plasmid['hits']])}\n")

            if(plasmid_id not in plasmid_dict):
                plasmid_dict[plasmid_id] = [0] * len(cfg.genome_path)
                log.info('Plasmid added: id=%s', plasmid_id)
            plasmid_dict[plasmid_id][genome_index] = 1

    with cfg.summary_path.open('w') as fh:
        fh.write(f'# {len(cfg.genome_path)} draft genome(s), {len(reference_plasmids)} reference plasmid(s)\n')
        fh.write('Genome\tPlasmid\tCoverage\tIdentity\tContigs\tContig IDs\n')
        for line in plasmid_string_summary:
            fh.write(line)

    if(plasmid_dict):
        write_cohort_table(plasmid_dict)

    if(plasmid_detected):
        json_path = cfg.output_path.joinpath("plasmids.json")
        tio.export_json(plasmid_detected, json_path)

def pooling(genome, reference_plasmids, index, lookup_plasmid_ids):
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
    filtered_hits = tb.filter_contig_hits(sample, hits, reference_plasmids, lookup_plasmid_ids)  # plasmid hits filtered by coverage and identity
    detected_plasmids = tp.detect_reference_plasmids(sample, filtered_hits, reference_plasmids)  # detect reference plasmids above cov/id thresholds

    # Write output files
    sample_summary_path = cfg.output_path.joinpath(f'{sample}-summary.tsv')
    with sample_summary_path.open('w') as ssp:
        ssp.write("plasmid\tcontig\tcontig start\tcontig end\tcontig length\tcoverage[%]\tidentity[%]\talignment length\tstrand\tplasmid start\tplasmid end\tplasmid length\n")

        for plasmid in detected_plasmids:

            plasmid_contigs_sorted = tp.reconstruct_plasmid(plasmid, contigs)

            prefix = f"{cfg.prefix}-{sample}-{plasmid['reference']}" if cfg.prefix else f"{sample}-{plasmid['reference']}"
            plasmid_contigs_path = cfg.output_path.joinpath(f'{prefix}-contigs.fna')
            plasmid_pseudosequence_path = cfg.output_path.joinpath(f'{prefix}-pseudo.fna')
            tio.export_sequences(plasmid_contigs_sorted, plasmid_contigs_path, description=True, wrap=True)
            tio.export_sequences([plasmid], plasmid_pseudosequence_path, description=True, wrap=True)

            # Write detailed plasmid hits to sample summary file
            for hit in plasmid['hits']:
                ssp.write(f"{plasmid['reference']}\t{hit['contig_id']}\t{hit['contig_start']}\t{hit['contig_end']}\t{hit['contig_length']}\t{hit['coverage']:.3f}\t{hit['perc_identity']:.3f}\t{hit['length']}\t{hit['strand']}\t{hit['reference_plasmid_start']}\t{hit['reference_plasmid_end']}\t{plasmid['length']}\n")

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

    return index, detected_plasmids


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
