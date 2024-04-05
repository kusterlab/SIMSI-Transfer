import logging
import multiprocessing
from pathlib import Path
import argparse

import pandas as pd

from . import __version__, __copyright__

logger = logging.getLogger(__name__)


class ArgumentParserWithLogger(argparse.ArgumentParser):
    def error(self, message):
        logger.error(f"Error parsing input arguments: {message}")
        super().error(message)


def parse_args(argv):
    desc = f'SIMSI-Transfer version {__version__}\n{__copyright__}' 
    apars = ArgumentParserWithLogger(
        description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    apars.add_argument('--mq_txt_folder', default=None, metavar="DIR",
                       help='''
                       Path to MaxQuant combined/txt output folder.
                       ''')

    apars.add_argument('--raw_folder', default=None, metavar="DIR",
                       help='''
                       Full path to folder containing .raw or .mzML files.
                       ''')

    apars.add_argument('--meta_input_file', default=None, metavar="DIR",
                       help='''
                       Tab separated file with a header line followed by rows containing mq_txt_folder, raw_folder and, optionally, tmt_reporter_correction_file.
                       ''')

    apars.add_argument('--stringencies', default="20,15,10", metavar="S",
                       help='''
                       Clustering thresholds at which to produce cluster files, listed as comma separated list. 
                       The higher the stringency value, the more strict the clustering is.
                       ''')

    apars.add_argument('--output_folder', type=Path, default=Path("./simsi_output"), metavar="DIR",
                       help='''
                       Full path to SIMSI output folder. This folder will contain the 
                       MaRaCluster results and a "summaries" folder which contains the 
                       SIMSI results in MaxQuant output format.
                       ''')
    
    apars.add_argument('--cache_folder', type=Path, default=Path("./simsi_output"), metavar="DIR",
                       help='''
                       Full path to SIMSI cache folder. The cache folder stores the mzML, 
                       extracted TMT intensities and MaRaCluster dat files such that 
                       they can be reused in multiple SIMSI runs.
                       ''')

    apars.add_argument('--num_threads', type=int, default=min(multiprocessing.cpu_count(), 4), metavar='N',
                       help='''
                       Number of threads, by default this is equal to the number of CPU cores available on the device.
                       ''')

    apars.add_argument('--tmt_reporter_correction_file', default="", metavar="DIR",
                       help='''
                       Path to TMT correction factor file, as exported from MaxQuant.
                       ''')

    apars.add_argument('--tmt_ms_level', default="ms2", metavar="S",
                       help='''
                       MS level of TMT quantification, either "ms2" or "ms3"
                       ''')

    apars.add_argument('--tmt_requantify', default=False, action='store_true',
                       help='''
                       Re-quantifies the TMT reporter ions directly from the raw file. 
                       This solves problems where MaxQuant assigns TMT reporter ions to the wrong scan in msmsScans.txt.
                       This appears especially problematic in MS3 data.
                       ''')

    apars.add_argument('--filter_decoys', default=False, action='store_true',
                       help='''
                       Removes decoys from MaxQuant results before PSM transfer.
                       ''')

    apars.add_argument('--skip_annotated_clusters', default=False, action='store_true',
                       help='''
                       Skip writing the p<x>_annotated_clusters.txt output files.
                       ''')
    
    apars.add_argument('--skip_msmsscans', default=False, action='store_true',
                       help='''
                       Skip writing the p<x>_msmsScans.txt output files.
                       ''')
    
    apars.add_argument('--skip_msms', default=False, action='store_true',
                       help='''
                       Skip writing the p<x>_msms.txt output files.
                       ''')
    
    apars.add_argument('--skip_evidence', default=False, action='store_true',
                       help='''
                       Skip writing the p<x>_evidence.txt output files.
                       ''')
    
    apars.add_argument('--ambiguity_decision', default='majority', metavar="S",
                       help='''
                       Determines which modified sequence to set for PTM-ambiguous clusters: 
                       'majority' returns the most frequently observed modified sequence in the cluster.
                       'all' returns unmodified sequence with all observed PTM positions.
                       ''')

    apars.add_argument('--add_plotting_columns', default=False, action='store_true',
                       help='''
                       Retains columns that might be needed for further processing in e.g. TMT curve plotting tools.
                       ''')

    apars.add_argument('--maximum_pep', type=int, default=100, metavar='N',
                       help='''
                       Maximum Posterior Error Probability (PEP) in percent of PSMs to be considered for transfers.
                       ''')

    # ------------------------------------------------
    args = apars.parse_args(argv)

    return args


def get_input_folders(args):
    if args.meta_input_file:
        if args.raw_folder:
            raise ValueError("Cannot use the --raw_folder and --meta_input_file parameters at the same time.")
        if args.mq_txt_folder:
            raise ValueError("Cannot use the --mq_txt_folder and --meta_input_file parameters at the same time.")

        meta_input_df = pd.read_csv(args.meta_input_file, sep='\t')
        meta_input_df.columns = ['mq_txt_folder', 'raw_folder', 'tmt_correction_file'][:len(meta_input_df.columns)]
        if 'tmt_correction_file' not in meta_input_df.columns:
            meta_input_df['tmt_correction_file'] = args.tmt_reporter_correction_file
    else:
        if not args.raw_folder:
            raise ValueError("Missing --raw_folder argument")
        if not args.mq_txt_folder:
            raise ValueError("Missing --mq_txt_folder argument")

        meta_input_df = pd.DataFrame(
            list(zip([args.mq_txt_folder], [args.raw_folder], [args.tmt_reporter_correction_file])),
            columns=['mq_txt_folder', 'raw_folder', 'tmt_correction_file'])

    return meta_input_df


def parse_stringencies(stringencies):
    if stringencies == '':
        stringencies = [20, 15, 10]
    else:
        try:
            stringencies = [int(i) for i in stringencies.split(',')]
        except ValueError:
            raise ValueError(
                'This is not a stringency list. Please input the list in the following format: 30,25,20,15,10,5')

    return stringencies


def parse_tmt_ms_level(extraction_level: str):
    if extraction_level == "ms2":
        return 2
    elif extraction_level == "ms3":
        return 3
    else:
        raise ValueError("--tmt_ms_level should be either ms2 or ms3.")


if __name__ == '__main__':
    raise NotImplementedError('Do not run this script.')
