import logging
import sys
import concurrent.futures as cf
import numpy as np

import tadrep.config as cfg
import tadrep.io as tio
import tadrep.blast as tb
import tadrep.plasmids as tp


log = logging.getLogger('DETECTION')


def detect_and_reconstruct():
    ############################################################################
    # Import plasmid sequences
    # - write multi Fasta file
    ############################################################################

    if(not cfg.db_data.get('plasmids', {})):
        log.debug("No plasmids in %s !", cfg.db_path)
        sys.exit(f"ERROR: No plasmids in database {cfg.db_path}!")

    if(not cfg.db_data.get('clusters', [])):
        log.debug("No Clusters in %s!", cfg.db_path)
        sys.exit(f"ERROR: No cluster in database {cfg.db_path}")

    log.info("Loaded %d cluster with %d plasmids", len(cfg.db_data['clusters']), len(cfg.db_data['plasmids'].keys()))
    cfg.verbose_print("Loaded data:")
    cfg.verbose_print(f"\t{len(cfg.db_data['clusters'])} cluster")
    cfg.verbose_print(f"\t{len(cfg.db_data['plasmids'].keys())} plasmids total")

    # Merge clusters and representative info from DB
    #reference_plasmids = {cluster['id']: cfg.db_data['plasmids'][cluster['representative']] for cluster in cfg.db_data['clusters']}
    reference_plasmids = {}
    for cluster in cfg.db_data['clusters']:
        reference_plasmid = cfg.db_data['plasmids'][cluster['representative']].copy()
        reference_plasmid.update(cluster)
        reference_plasmids[cluster['id']] = reference_plasmid
    references_path = cfg.output_path.joinpath('references.fna')
    tio.export_sequences(reference_plasmids.values(), references_path)

    cfg.verbose_print(f"Found {len(reference_plasmids)} representative plasmid(s)")
    log.info("Found %d representative plasmid(s)", len(reference_plasmids))



    ############################################################################
    # Prepare summary output file
    # - create file
    # - write header into file
    ############################################################################
    plasmid_dict = {}
    plasmid_string_summary = []
    plasmids_detected = {}

    cfg.verbose_print('Analyze genome sequences...')
    futures = []
    with cf.ThreadPoolExecutor(max_workers=cfg.threads) as pool:
        for index, genome_path in enumerate(cfg.genome_path):
            futures.append(pool.submit(detect_plasmids, genome_path, reference_plasmids, index))

    for future in futures:
        genome_index, plasmid_summary = future.result()
        for plasmid in plasmid_summary:
            reference_id = plasmid['reference']
            if(reference_id not in plasmids_detected):
                plasmids_detected[reference_id] = {k: plasmid[k] for k in ['id', 'reference', 'length']}
                plasmids_detected[reference_id]['found_in'] = {}
            plasmids_detected[reference_id]['found_in'][plasmid['genome']] = plasmid['hits']

            # Create string for plasmid summary
            plasmid_string_summary.append(f"{plasmid['genome']}\t{plasmid['reference']}\t{plasmid['coverage']:.3f}\t{plasmid['identity']:.3f}\t{len(plasmid['hits'])}\t{','.join([hit['contig_id'] for hit in plasmid['hits']])}\n")

            if(reference_id not in plasmid_dict):
                plasmid_dict[reference_id] = [0] * len(cfg.genome_path)
                log.info('Plasmid added: id=%s', reference_id)
            plasmid_dict[reference_id][genome_index] = 1

    with cfg.summary_path.open('w') as fh:
        fh.write(f'# {len(cfg.genome_path)} draft genome(s), {len(reference_plasmids)} reference plasmid(s)\n')
        fh.write('Genome\tPlasmid\tCoverage\tIdentity\tContigs\tContig IDs\n')
        for line in plasmid_string_summary:
            fh.write(line)

    if(plasmid_dict):
        write_cohort_table(plasmid_dict)
        write_plasmids_info(plasmid_dict, reference_plasmids)

    if(plasmids_detected):
        for reference_id, plasmid_data in plasmids_detected.items():
            for cluster in cfg.db_data['clusters']:
                if cluster['id'] == reference_id:
                    cluster['found_in'] = plasmid_data['found_in']
        tio.export_json(cfg.db_data, cfg.db_path)


def detect_plasmids(genome, reference_plasmids, index):
    log_pool = logging.getLogger('PROCESS')

    # Import draft genome contigs
    try:
        contigs = tio.import_sequences(genome, sequence=True)
        log_pool.info('imported genome contigs: genome=%s, # contigs=%i', genome, len(contigs))
    except ValueError:
        log_pool.error('wrong genome file format!', exc_info=True)
        sys.exit('ERROR: wrong genome file format!')

    sample = genome.stem
    blast_output_path = cfg.tmp_path.joinpath(f'{sample}-{index}-blastn.tsv')
    hits = tb.search_contigs(genome, blast_output_path)  # plasmid raw hits
    filtered_hits = tb.filter_contig_hits(sample, hits, reference_plasmids)  # plasmid hits filtered by coverage and identity
    detected_plasmids = tp.detect_reference_plasmids(sample, filtered_hits, reference_plasmids)  # detect reference plasmids above cov/id thresholds

    # Write output files
    sample_summary_path = cfg.output_path.joinpath(f'{sample}-summary.tsv')
    with sample_summary_path.open('w') as ssp:
        ssp.write("plasmid\tcontig\tcontig start\tcontig end\tcontig length\tcoverage[%]\tidentity[%]\talignment length\tstrand\tplasmid start\tplasmid end\tplasmid length\n")

        for plasmid in detected_plasmids:
            plasmid_contigs_sorted = tp.reconstruct_plasmid(plasmid, contigs)
            prefix = f"{cfg.prefix}-{sample}-{plasmid['reference']}" if cfg.prefix else f"{sample}-{plasmid['reference']}"
            plasmid_contigs_path = cfg.output_path.joinpath(f'{prefix}-contigs.fna')
            tio.export_sequences(plasmid_contigs_sorted, plasmid_contigs_path, description=True, wrap=True)
            plasmid_pseudosequence_path = cfg.output_path.joinpath(f'{prefix}-pseudo.fna')
            tio.export_sequences([plasmid], plasmid_pseudosequence_path, description=True, wrap=True)

            # Write detailed plasmid hits to sample summary file
            for hit in plasmid['hits']:
                ssp.write(f"{plasmid['reference']}\t{hit['contig_id']}\t{hit['contig_start']}\t{hit['contig_end']}\t{hit['contig_length']}\t{hit['coverage']:.3f}\t{hit['perc_identity']:.3f}\t{hit['length']}\t{hit['strand']}\t{hit['reference_plasmid_start']}\t{hit['reference_plasmid_end']}\t{plasmid['length']}\n")

    cfg.lock.acquire()
    log_pool.debug('lock acquired: genome=%s, index=%s', sample, index)
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
    log_pool.debug('lock released: genome=%s, index=%s', sample, index)
    cfg.lock.release()

    return index, detected_plasmids


def write_cohort_table(plasmid_dict):
    plasmid_cohort_path = cfg.output_path.joinpath('plasmids.distribution.tsv')
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
                fh.write(f'\t{"1" if plasmid else "0"}')
            fh.write('\n')


def write_plasmids_info(plasmid_dict, reference_plasmids):
    plasmid_info_path = cfg.output_path.joinpath('plasmids.info.tsv')
    with plasmid_info_path.open('w') as fh:
        fh.write('Plasmid\tRepresentative\tLength\tGC\tCDS\tINC_Types\n')
        for plasmid_id in plasmid_dict.keys():
            plasmid_inc_types = ', '.join([inc_type['type'] for inc_type in reference_plasmids[plasmid_id]['inc_types']]) if len(reference_plasmids[plasmid_id]["inc_types"]) > 0 else "-"
            fh.write(f"{reference_plasmids[plasmid_id]['id']}\t{reference_plasmids[plasmid_id]['representative']}\t{reference_plasmids[plasmid_id]['length']}\t{reference_plasmids[plasmid_id]['gc_content']}\t{len(reference_plasmids[plasmid_id]['cds'])}\t{plasmid_inc_types:>}\n")
