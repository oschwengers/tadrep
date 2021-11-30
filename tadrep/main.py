import logging
import sys
import shutil
import numpy as np
from pathlib import Path
import multiprocessing as mp

import tadrep
import tadrep.blast as tb
import tadrep.config as cfg
import tadrep.fasta as fasta
import tadrep.plasmids as tp
import tadrep.utils as tu
import tadrep.visuals as tv


def main():
    # parse arguments
    args = tu.parse_arguments()

    ############################################################################
    # Setup logging
    ############################################################################
    try:
        output_path = Path(args.output) if args.output else Path.cwd()
        if(not output_path.exists()):
            output_path.mkdir(parents=True, exist_ok=True)
        output_path = output_path.resolve()
        cfg.output_path = output_path
    except:
        sys.exit(f'ERROR: could not resolve or create output directory ({args.output})!')
    log_prefix = args.prefix if args.prefix else "tadrep"
    logging.basicConfig(
        filename=str(output_path.joinpath(f'{log_prefix}.log')),
        filemode='w',
        format='%(asctime)s - %(processName)s - %(levelname)s - %(name)s - %(message)s',
        level=logging.DEBUG if args.verbose else logging.INFO
    )
    log = logging.getLogger('MAIN')
    log.info('version %s', tadrep.__version__)
    log.info('command line: %s', ' '.join(sys.argv))


    ############################################################################
    # Checks and configurations
    # - check parameters and setup global configuration
    # - test database
    # - test binary dependencies
    ############################################################################
    cfg.setup(args)  # check parameters and prepare global configuration

    verboseprint = print if cfg.verbose else lambda *a, **k: None
    verboseprint(f'TaDReP v{tadrep.__version__}')
    verboseprint('Options and arguments:')
    verboseprint(f'\tinput: {[str(genome) for genome in cfg.genome_path]}')
    verboseprint(f'\tplasmid(s): {cfg.plasmids_path}')
    verboseprint(f'\tdatabase path: {cfg.database_path}')
    verboseprint(f'\toutput: {cfg.output_path}')
    verboseprint(f'\tprefix: {cfg.prefix}')
    verboseprint(f'\ttmp directory: {cfg.tmp_path}')
    verboseprint(f'\t# threads: {cfg.threads}')


    ############################################################################
    # Import plasmid sequences
    # - parse contigs in Fasta file
    ############################################################################
    verboseprint('parse plasmids sequences...')
    if(cfg.plasmids_path):
        try:
            plasmids = fasta.import_sequences(cfg.plasmids_path)
            log.info('imported plasmids: sequences=%i, file=%s', len(plasmids), cfg.plasmids_path)
            verboseprint(f'\timported: {len(plasmids)}')
        except ValueError:
            log.error('wrong plasmids file format!', exc_info=True)
            sys.exit('ERROR: wrong plasmids file format!')
    else:
        try:
            plasmids = fasta.import_tsv(cfg.database_path)
            log.info('imported plasmids: sequence=%i, file=%s', len(plasmids), cfg.database_path)
            verboseprint(f'\timported: {len(plasmids)}')
        except:
            log.error('wrong database path!', exc_info=True)
            sys.exit('ERROR: wrong database path!')
    ############################################################################
    # Prepare summary output file
    # - create file
    # - write header into file
    ############################################################################
    plasmid_list = {}
    plasmid_string_summary = []

    values = ((genome, plasmids, genome_index) for genome_index, genome in enumerate(cfg.genome_path))
    with mp.Pool(cfg.threads) as pool:
        genomes_summary = pool.starmap(pooling, values)

    for genome in genomes_summary:
        genome_index = genome[0]
        plasmid_summary = genome[1]
        for plasmid_id, summary_string in plasmid_summary.items():
            plasmid_string_summary.append(summary_string)
            if(not plasmid_list.get(plasmid_id, None)):
                plasmid_list[plasmid_id] = [0] * len(cfg.genome_path)
                log.info('Plasmid added: id=%s', plasmid_id)
            plasmid_list[plasmid_id][genome_index] = 1

    cfg.summary_path = cfg.output_path.joinpath('summary.tsv')
    with cfg.summary_path.open('w') as fh:
        fh.write(f'# {len(cfg.genome_path)} draft genome(s), {len(plasmids)} plasmid(s)\n')
        fh.write('file\tplasmid_id\tcoverage\tidentity\tcontigs\tcontig_ids\n')
        for line in plasmid_string_summary:
            fh.write(line)

    if(plasmid_list):
        write_cohort_table(plasmid_list)

    # remove tmp dir
    shutil.rmtree(str(cfg.tmp_path))
    log.debug('removed tmp dir: %s', cfg.tmp_path)


def pooling(genome, plasmids, index):
    log = logging.getLogger('Process')
    verboseprint = print if cfg.verbose else lambda *a, **k: None

    verboseprint('parse genome sequences...')
    # Import draft genome contigs
    try:
        contigs = fasta.import_sequences(genome, sequence=True)
        log.info('imported genomes: sequences=%i, file=%s', len(contigs), genome)
        verboseprint(f'\timported: {len(contigs)}, file: {genome.stem}')
    except ValueError:
        log.error('wrong genome file format!', exc_info=True)
        sys.exit('ERROR: wrong genome file format!')

    sample = genome.stem
    blast_output_path = cfg.tmp_path.joinpath(f'{sample}-{index}-blastn.tsv')
    log.debug('Blast output: path=%s', blast_output_path)

    verboseprint('Blasting against reference...')
    hits = tb.search_contigs(genome, blast_output_path)  # List of BLAST hits
    verboseprint('Filtering hits...')
    filtered_hits = tb.filter_contig_hits(hits)  # Dictionary of Plasmids, with Blast hits filtered by coverage and identity
    verboseprint('Detecting plasmids...')
    detected_plasmids = tp.detect_plasmids(filtered_hits, plasmids)  # List of detect reference plasmids above cov/id thresholds

    # Write output files
    verboseprint('write genome sequences...')

    sample_summary_path = cfg.output_path.joinpath(f'{sample}-summary.tsv')
    cfg.lock.acquire()
    log.debug('Lock acquired: genome=%s, index=%s', sample, index)
    if(sample_summary_path.is_file()):
        sample_summary_path.unlink()
    log.info('Summary file: sample=%s, path=%s', sample, sample_summary_path)
    log.debug('Lock released: genome=%s, index=%s', sample, index)
    cfg.lock.release()

    plasmid_summary_strings = {}

    with sample_summary_path.open('w+') as ssp:
        ssp.write(f"plasmid\tcontig\tcontig start\tcontig end\tcontig length\tcoverage[%]\tidentity[%]\talignment length\tstrand\tplasmid start\tplasmid end\tplasmid length\n")

        for plasmid in detected_plasmids:

            prefix = f"{cfg.prefix}-{sample}-{plasmid['id']}" if cfg.prefix else f"{sample}-{plasmid['id']}"
            plasmid_contigs_path = cfg.output_path.joinpath(f'{prefix}-contigs.fna')
            plasmid_pseudosequence_path = cfg.output_path.joinpath(f'{prefix}-pseudo.fna')

            log.debug('prepare output: plasmid-id=%s, contigs-path=%s, assembly-path=%s', plasmid['id'], plasmid_contigs_path, plasmid_pseudosequence_path)
            plasmid_contigs_sorted = tp.reconstruct_plasmid(plasmid, sample, contigs)

            fasta.export_sequences(plasmid_contigs_sorted, plasmid_contigs_path, description=True, wrap=True)
            fasta.export_sequences([plasmid], plasmid_pseudosequence_path, description=True, wrap=True)

            # Create string for plasmid summary
            plasmid_string = [
                f'{sample}',
                f"{plasmid['id']}",
                f"{plasmid['coverage']:.3f}",
                f"{plasmid['identity']:.3f}",
                f"{len(plasmid['hits'])}",
                f"{','.join([hit['contig_id'] for hit in plasmid['hits']])}\n"
            ]
            plasmid_summary_strings[plasmid['id']] = '\t'.join(plasmid_string)

            # Write detailed plasmid hits to sample summary file
            for hit in plasmid['hits']:
                ssp.write(f"{plasmid['id']}\t")
                ssp.write(f"{hit['contig_id']}\t")
                ssp.write(f"{hit['contig_start']}\t")
                ssp.write(f"{hit['contig_end']}\t")
                ssp.write(f"{hit['contig_length']}\t")
                ssp.write(f"{hit['coverage']:.3f}\t")
                ssp.write(f"{hit['perc_identity']:.3f}\t")
                ssp.write(f"{hit['length']}\t")
                ssp.write(f"{hit['strand']}\t")
                ssp.write(f"{hit['plasmid_start']}\t")
                ssp.write(f"{hit['plasmid_end']}\t")
                ssp.write(f"{plasmid['length']}\n")

            png_path = cfg.output_path.joinpath(f'{sample}-{plasmid["id"]}.png')
            tv.create_plasmid_figure(plasmid, genome.name, png_path)

    cfg.lock.acquire()
    log.debug('Lock acquired: genome=%s, index=%s', sample, index)
    print(f'\nGenome file: {sample}, total contigs: {len(contigs)}, found plasmids: {len(detected_plasmids)}\n')
    for plasmid in detected_plasmids:
        # Create console output
        if (cfg.verbose):
            print(f"\tplasmid: {plasmid['id']}\t({plasmid['length']} bp)\tcontig hits = {len(plasmid['hits'])}\t coverage = {plasmid['coverage'] * 100:1.1f}%\tidentity = {plasmid['identity'] * 100:1.1f}%")
            print(f"\t{'contig_id':^17} hit_length contig_length contig_start contig_end strand plasmid_start plasmid_end coverage[%] identity[%]")
            for hit in plasmid['hits']:
                contig = contigs[hit['contig_id']]
                print(f"\t{hit['contig_id']:^17} {hit['length']:>10} {contig['length']:>13} {hit['contig_start']:>12} {hit['contig_end']:>10} {hit['strand']:^6} {hit['plasmid_start']:>13} {hit['plasmid_end']:>11} {(hit['length'] / contig['length'] * 100):>8.1f}    {hit['perc_identity'] * 100:>8.1f}")
            print('\n')
        else:
            print(f"{sample}\t{plasmid['id']}\t{plasmid['length']}\t{len(plasmid['hits'])}\t{plasmid['coverage']:f}\t{plasmid['identity']:f}\t{','.join(contig['contig_id'] for contig in plasmid['hits'])}")
    log.debug('Lock released: genome=%s, index=%s', sample, index)
    cfg.lock.release()

    return index, plasmid_summary_strings


def write_cohort_table(plasmid_list):
    # Write cohort table
    plasmid_cohort_path = cfg.output_path.joinpath('plasmids.tsv')
    if(plasmid_cohort_path.is_file()):
        plasmid_cohort_path.unlink()
    with plasmid_cohort_path.open('w+') as fh:
        # Write plasmid header
        plasmid_order = []
        for plasmid in plasmid_list.keys():
            plasmid_order.append(plasmid_list[plasmid])
            fh.write(f'\t{plasmid}')
        fh.write('\n')
        # Transpose information for easier writing
        transposed_plasmid_order = np.array(plasmid_order).T.tolist()
        # Mark which plasmid was found for each draft genome
        for num_genome, genome in enumerate(cfg.genome_path):
            sample = genome.stem
            fh.write(f'{sample}')
            for plasmid in transposed_plasmid_order[num_genome]:
                fh.write(f'\t{"X" if plasmid else "-"}')
            fh.write('\n')


if __name__ == '__main__':
    main()
