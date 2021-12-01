import logging
import multiprocessing as mp
import numpy as np
import shutil
import sys

from pathlib import Path

import tadrep
import tadrep.blast as tb
import tadrep.config as cfg
import tadrep.io as tio
import tadrep.plasmids as tp
import tadrep.utils as tu
import tadrep.visuals as tv


def main():
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
    verboseprint(f"\tgenome(s): {', '.join([genome.name for genome in cfg.genome_path])}")
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
    if(cfg.plasmids_path):
        try:
            verboseprint('\nparse plasmids sequences...')
            plasmids = tio.import_sequences(cfg.plasmids_path)
            log.info('imported plasmids: sequences=%i, file=%s', len(plasmids), cfg.plasmids_path)
            verboseprint(f'\timported: {len(plasmids)}')
        except ValueError:
            log.error('wrong plasmids file format!', exc_info=True)
            sys.exit('ERROR: wrong plasmids file format!')
    else:
        try:
            verboseprint('\nload plasmid database...')
            plasmids = tio.import_tsv(cfg.database_path)
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
    plasmid_dict = {}
    plasmid_string_summary = []

    verboseprint('Analyze genome sequences...')
    values = ((genome_path, plasmids, genome_index) for genome_index, genome_path in enumerate(cfg.genome_path))
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
        fh.write(f'# {len(cfg.genome_path)} draft genome(s), {len(plasmids)} plasmid(s)\n')
        fh.write('Genome\tPlasmid\tCoverage\tIdentity\tContigs\tContig IDs\n')
        for line in plasmid_string_summary:
            fh.write(line)

    if(plasmid_dict):
        write_cohort_table(plasmid_dict)

    # remove tmp dir
    shutil.rmtree(str(cfg.tmp_path))
    log.debug('removed tmp dir: %s', cfg.tmp_path)


def pooling(genome, plasmids, index):
    log = logging.getLogger('Process')

    # Import draft genome contigs
    try:
        contigs = tio.import_sequences(genome, sequence=True)
        log.info('imported genomes: sequences=%i, file=%s', len(contigs), genome)
    except ValueError:
        log.error('wrong genome file format!', exc_info=True)
        sys.exit('ERROR: wrong genome file format!')

    sample = genome.stem
    blast_output_path = cfg.tmp_path.joinpath(f'{sample}-{index}-blastn.tsv')
    log.debug('Blast output: path=%s', blast_output_path)

    hits = tb.search_contigs(genome, blast_output_path)  # List of BLAST hits
    filtered_hits = tb.filter_contig_hits(hits)  # Dictionary of Plasmids, with Blast hits filtered by coverage and identity
    detected_plasmids = tp.detect_plasmids(filtered_hits, plasmids)  # List of detect reference plasmids above cov/id thresholds

    # Write output files
    sample_summary_path = cfg.output_path.joinpath(f'{sample}-summary.tsv')
    plasmid_summary_strings = {}
    with sample_summary_path.open('w') as ssp:
        ssp.write(f"plasmid\tcontig\tcontig start\tcontig end\tcontig length\tcoverage[%]\tidentity[%]\talignment length\tstrand\tplasmid start\tplasmid end\tplasmid length\n")

        for plasmid in detected_plasmids:

            prefix = f"{cfg.prefix}-{sample}-{plasmid['reference']}" if cfg.prefix else f"{sample}-{plasmid['reference']}"
            plasmid_contigs_path = cfg.output_path.joinpath(f'{prefix}-contigs.fna')
            plasmid_pseudosequence_path = cfg.output_path.joinpath(f'{prefix}-pseudo.fna')

            log.debug('prepare output: plasmid-id=%s, contigs-path=%s, assembly-path=%s', plasmid['id'], plasmid_contigs_path, plasmid_pseudosequence_path)
            plasmid_contigs_sorted = tp.reconstruct_plasmid(plasmid, sample, contigs)

            tio.export_sequences(plasmid_contigs_sorted, plasmid_contigs_path, description=True, wrap=True)
            tio.export_sequences([plasmid], plasmid_pseudosequence_path, description=True, wrap=True)

            # Create string for plasmid summary
            plasmid_summary_strings[plasmid['reference']] = f"{sample}\t{plasmid['reference']}\t{plasmid['coverage']:.3f}\t{plasmid['identity']:.3f}\t{len(plasmid['hits'])}\t{','.join([hit['contig_id'] for hit in plasmid['hits']])}\n"

            # Write detailed plasmid hits to sample summary file
            for hit in plasmid['hits']:
                ssp.write(f"{plasmid['reference']}\t{hit['contig_id']}\t{hit['contig_start']}\t{hit['contig_end']}\t{hit['contig_length']}\t{hit['coverage']:.3f}\t{hit['perc_identity']:.3f}\t{hit['length']}\t{hit['strand']}\t{hit['plasmid_start']}\t{hit['plasmid_end']}\t{plasmid['length']}\n")

            png_path = cfg.output_path.joinpath(f"{sample}-{plasmid['reference']}.png")
            tv.create_plasmid_figure(plasmid, genome.name, png_path)

    cfg.lock.acquire()
    log.debug('Lock acquired: genome=%s, index=%s', sample, index)
    print(f'\n\nGenome: {sample}, contigs: {len(contigs)}, detected plasmids: {len(detected_plasmids)}')
    for plasmid in detected_plasmids:
        # Create console output
        if (cfg.verbose):
            print(f"\n\tplasmid: {plasmid['reference']}, length: {plasmid['length']} bp, contig hits: {len(plasmid['hits'])}, coverage: {plasmid['coverage'] * 100:1.1f}%, identity: {plasmid['identity'] * 100:1.1f}%")
            print(f"\t{'contig':^17}  alignment length  contig length  contig start  contig end  strand  plasmid start  plasmid end  coverage[%]  identity[%]")
            for hit in plasmid['hits']:
                contig = contigs[hit['contig_id']]
                print(f"\t{hit['contig_id']:^18} {hit['length']:>11} {contig['length']:>14} {hit['contig_start']:>13} {hit['contig_end']:>13} {hit['strand']:^13} {hit['plasmid_start']:>9} {hit['plasmid_end']:>12} {(hit['length'] / contig['length'] * 100):>13.1f} {hit['perc_identity'] * 100:>12.1f}")
        else:
            print(f"{sample}\t{plasmid['id']}\t{plasmid['length']}\t{len(plasmid['hits'])}\t{plasmid['coverage']:f}\t{plasmid['identity']:f}\t{','.join(contig['contig_id'] for contig in plasmid['hits'])}")
    log.debug('Lock released: genome=%s, index=%s', sample, index)
    cfg.lock.release()

    return index, plasmid_summary_strings


def write_cohort_table(plasmid_dict):
    plasmid_cohort_path = cfg.output_path.joinpath('plasmids.tsv')
    with plasmid_cohort_path.open('w') as fh:
        # Write plasmid header
        plasmid_order = []
        for plasmid in plasmid_dict.keys():
            plasmid_order.append(plasmid_dict[plasmid])
            fh.write(f'\t{plasmid}')
        fh.write('\n')
        
        transposed_plasmid_order = np.array(plasmid_order).T.tolist()  # Transpose information for easier writing
        for num_genome, genome in enumerate(cfg.genome_path):  # Mark which plasmid was found for each draft genome
            sample = genome.stem
            fh.write(f'{sample}')
            for plasmid in transposed_plasmid_order[num_genome]:
                fh.write(f'\t{"X" if plasmid else "-"}')
            fh.write('\n')


if __name__ == '__main__':
    main()
