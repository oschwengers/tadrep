import logging
import shutil
import sys

from pathlib import Path

import tadrep
import tadrep.config as cfg
import tadrep.utils as tu
import tadrep.setup as ts
import tadrep.extract as te
import tadrep.cluster as tcl
import tadrep.characterize as tc
import tadrep.detect as td
import tadrep.visualize as tv
import database.main as dm


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

    cfg.verbose_print(f'TaDReP v{tadrep.__version__}')
    cfg.verbose_print('Options and arguments:')
    cfg.verbose_print(f'\toutput: {cfg.output_path}')
    cfg.verbose_print(f'\tprefix: {cfg.prefix}')
    cfg.verbose_print(f'\ttmp directory: {cfg.tmp_path}')
    cfg.verbose_print(f'\t# threads: {cfg.threads}')

    if(args.subcommand == 'database'):
        cfg.setup_database(args)
        dm.create_database()

    if(args.subcommand == 'setup'):
        print('\nSetup started...')
        ts.download_inc_reference()

    if(args.subcommand == "extract"):
        print('\nExtraction started...')
        cfg.setup_extract(args)
        te.extract()

    elif(args.subcommand == "characterize"):
        print('\nCharacterization started...')
        cfg.setup_characterize(args)
        tc.characterize()

    elif(args.subcommand == "cluster"):
        print('\nClustering started...')
        cfg.setup_cluster(args)
        tcl.cluster_plasmids()

    elif(args.subcommand == "detect"):
        cfg.setup_detect(args)
        print(f"\tgenome(s): {', '.join([genome.name for genome in cfg.genome_path])}")

        print('\nDetection and reconstruction started ...')
        td.detect_and_reconstruct()

    elif(args.subcommand == "visualize"):
        print('\nVisualization started...')
        cfg.setup_visualize(args)
        tv.plot()

    # remove tmp dir
    shutil.rmtree(str(cfg.tmp_path))
    log.debug('removed tmp dir: %s', cfg.tmp_path)


if __name__ == '__main__':
    main()
