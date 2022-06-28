import logging
import shutil
import sys

from pathlib import Path

import tadrep
import tadrep.config as cfg
import tadrep.utils as tu
import tadrep.detection as td
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
    log_prefix = args.prefix if args.prefix else 'tadrep'
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

    if(cfg.extraction):
        verboseprint('\nExtraction started...')

    if(cfg.characterization):
        verboseprint('\nCharacterization started...')

    if(cfg.clustering):
        verboseprint('\nClustering started...')

    if(cfg.detection):
        verboseprint('\nDetection and reconstruction started ...')
        td.detect_and_reconstruct()

    if(cfg.visualization):
        verboseprint('\nVisualization started...')
        tv.create_plots()

    # remove tmp dir
    shutil.rmtree(str(cfg.tmp_path))
    log.debug('removed tmp dir: %s', cfg.tmp_path)


if __name__ == '__main__':
    main()
