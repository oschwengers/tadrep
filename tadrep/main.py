import logging
import sys
import shutil
from pathlib import Path

import tadrep
import tadrep.blast as tb
import tadrep.config as cfg
import tadrep.fasta as fasta
import tadrep.plasmids as tp
import tadrep.utils as tu


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
        format='%(asctime)s - %(levelname)s - %(name)s - %(message)s',
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
    if(cfg.verbose):
        print(f'TaDReP v{tadrep.__version__}')
        print('Options and arguments:')
        print(f'\tinput: {cfg.genome_path}')
        print(f'\tplasmid(s): {cfg.plasmids_path}')
        print(f'\toutput: {cfg.output_path}')
        print(f'\tprefix: {cfg.prefix}')
        print(f'\ttmp directory: {cfg.tmp_path}')
        print(f'\t# threads: {cfg.threads}')

    ############################################################################
    # Import plasmid sequences
    # - parse contigs in Fasta file
    ############################################################################
    print('parse plasmids sequences...')
    try:
        plasmids = fasta.import_sequences(cfg.plasmids_path)
        log.info('imported plasmids: sequences=%i, file=%s', len(plasmids), cfg.plasmids_path)
        print(f'\timported: {len(plasmids)}')
    except ValueError:
        log.error('wrong plasmids file format!', exc_info=True)
        sys.exit('ERROR: wrong plasmids file format!')

    ############################################################################
    # Prepare summary output file
    # - create file
    # - write header into file
    ############################################################################
    cfg.summary_path = cfg.output_path.joinpath(f'{log_prefix}-summary.tsv')
    with cfg.summary_path.open('w') as fh:
        fh.write(f'# {len(cfg.genome_path)} draft genome(s), {len(plasmids)} plasmid(s)\n')
        fh.write('file\tplasmid_id\tcontig_id\tcontig_start\tcontig_end\tcontig_length\tstrand\tplasmid_start\tplasmid_end\tlength\tcoverage\tidentity\tevalue\n')

        ############################################################################
        # Import draft genome contigs
        # - parse contigs in Fasta file
        # - apply contig length filter
        ############################################################################
        print('parse genome sequences...')
        for genome in cfg.genome_path:
            try:
                contigs = fasta.import_sequences(genome)
                log.info('imported genomes: sequences=%i, file=%s', len(contigs), genome)
                print(f'\timported: {len(contigs)}, file: {genome.stem}')
            except ValueError:
                log.error('wrong genome file format!', exc_info=True)
                sys.exit('ERROR: wrong genome file format!')

            unfiltered_hits = tb.search_contigs(genome, cfg.plasmids_path)
            filtered_hits = tb.filter_contig_hits(unfiltered_hits)
            detected_plasmids = tp.detect_plasmids(filtered_hits, plasmids)

            ############################################################################
            # Write output files
            # - write comprehensive annotation results as JSON
            # - write optional output files in GFF3/GenBank/EMBL formats
            # - remove temp directory
            ############################################################################
            print('write genome sequences...\n')
            genome_file = genome.stem
            print(f'\n\nGenome file: {genome_file}, total contigs: {len(contigs)}, found plasmids: {len(detected_plasmids)}\n')
            for plasmid in detected_plasmids:
                prefix = f"{cfg.prefix}-{genome_file}-{plasmid['id']}" if cfg.prefix else f"{genome_file}-{plasmid['id']}"
                plasmid_contigs_path = cfg.output_path.joinpath(f'{prefix}-contigs.fna')
                plasmid_pseudosequence_path = cfg.output_path.joinpath(f'{prefix}-pseudo.fna')

                log.debug('prepare output: plasmid-id=%s, contigs-path=%s, assembly-path=%s', plasmid['id'], plasmid_contigs_path, plasmid_pseudosequence_path)
                matched_contigs_sorted = tp.reconstruct_plasmid(plasmid, genome_file, contigs)

                fasta.export_sequences(matched_contigs_sorted, plasmid_contigs_path, description=True, wrap=True)
                fasta.export_sequences([plasmid], plasmid_pseudosequence_path, description=True, wrap=True)

                reference_plasmid = plasmids[plasmid['id']]

                if (cfg.verbose):
                    print(f"\tplasmid: {plasmid['id']}\t({reference_plasmid['length']} bp)\tcontig hits = {len(plasmid['hits'])}\t coverage = {plasmid['coverage'] * 100:1.1f}%\tidentity = {plasmid['identity'] * 100:1.1f}%")
                    print(f"\t{'contig_id':^17} hit_length contig_length contig_start contig_end strand plasmid_start plasmid_end coverage[%] identity[%]")
                else:
                    print(f"{cfg.prefix}\t{plasmid['id']}\t{reference_plasmid['length']}\t{len(plasmid['hits'])}\t{plasmid['coverage']:f}\t{plasmid['identity']:f}\t{','.join(contig['contig_id'] for contig in plasmid['hits'])}")

                for hit in plasmid['hits']:
                    fh.write(f'{genome}\t')
                    fh.write(f"{plasmid['id']}\t")
                    fh.write(f"{hit['contig_id']}\t")
                    fh.write(f"{hit['contig_start']}\t")
                    fh.write(f"{hit['contig_end']}\t")
                    fh.write(f"{hit['contig_length']}\t")
                    fh.write(f"{hit['strand']}\t")
                    fh.write(f"{hit['plasmid_start']}\t")
                    fh.write(f"{hit['plasmid_end']}\t")
                    fh.write(f"{hit['length']}\t")
                    fh.write(f"{hit['coverage']:.3f}\t")
                    fh.write(f"{hit['perc_identity']:.3f}\t")
                    fh.write(f"{hit['evalue']}\n")

                    if (cfg.verbose):
                        contig = contigs[hit['contig_id']]
                        print(f"\t{hit['contig_id']:^17} {hit['length']:>10} {contig['length']:>13} {hit['contig_start']:>12} {hit['contig_end']:>10} {hit['strand']:^6} {hit['plasmid_start']:>13} {hit['plasmid_end']:>11} {(hit['length'] / contig['length'] * 100):>8.1f}    {hit['perc_identity'] * 100:>8.1f}")
                print('\n')

    # remove tmp dir
    shutil.rmtree(str(cfg.tmp_path))
    log.debug('removed tmp dir: %s', cfg.tmp_path)


if __name__ == '__main__':
    main()
