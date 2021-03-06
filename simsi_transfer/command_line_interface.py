import logging
import multiprocessing
from pathlib import Path
import argparse

import pandas as pd

logger = logging.getLogger(__name__)


class ArgumentParserWithLogger(argparse.ArgumentParser):
    def error(self, message):
        logger.error(f"Error parsing input arguments: {message}")
        super().error(message)


def parse_args(argv):
    apars = ArgumentParserWithLogger(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    apars.add_argument('--mq_txt_folder', default=None, metavar="DIR",
                       help='''Path to MaxQuant combined/txt output folder.''')

    apars.add_argument('--raw_folder', default=None, metavar="DIR",
                       help='''Full path to folder containing .raw or .mzML files.''')

    apars.add_argument('--meta_input_file', default=None, metavar="DIR",
                       help='''Tab separated file with a header line followed by rows containing mq_txt_folder, raw_folder and, optionally, tmt_reporter_correction_file.''')

    apars.add_argument('--stringencies', default="20,15,10", metavar="S",
                       help='''Clustering thresholds at which to produce cluster files, listed as comma separated list. 
                               The higher the stringency value, the more strict the clustering is.
                          ''')

    apars.add_argument('--output_folder', default="./simsi_output", metavar="DIR",
                       help='''Full path to desired SIMSI output folder.''')

    apars.add_argument('--num_threads', type=int, default=min(multiprocessing.cpu_count(), 4), metavar='N',
                       help='''Number of threads, by default this is equal to the number of CPU cores available on the device.
                          ''')

    apars.add_argument('--tmt_reporter_correction_file', default="", metavar="DIR",
                       help='''(optional) Path to TMT correction factor file, as exported from MaxQuant.''')

    apars.add_argument('--tmt_ms_level', default="ms2", metavar="S",
                       help='''MS level of TMT quantification, either "ms2" or "ms3"''')

    apars.add_argument('--tmt_requantify', default=False, action='store_true',
                       help='''Re-quantifies the TMT reporter ions directly from the raw file. 
                               This solves problems where MaxQuant assigns TMT reporter ions to the wrong scan in msmsScans.txt.
                               This appears especially problematic in MS3 data.''')

    apars.add_argument('--ambiguity_decision', default=False, action='store_true',
                       help='''Enables the decision for one spectrum whenever PTM-isomeric clusters are encountered.''')
    # TODO: Change into different kinds of handling; str for decision

    # ------------------------------------------------
    args = apars.parse_args(argv)
    
    raw_folders, mq_txt_folders, tmt_correction_files = get_input_folders(args)
    
    return mq_txt_folders, raw_folders, parse_stringencies(args.stringencies), Path(
        args.output_folder), args.num_threads, tmt_correction_files, args.tmt_ms_level, \
        args.tmt_requantify, args.ambiguity_decision, args.meta_input_file


def get_input_folders(args):
    if args.meta_input_file:
        if args.raw_folder:
            logging.error("Cannot use the --raw_folder and --meta_input_file parameters at the same time.")
        if args.mq_txt_folder:
            logging.error("Cannot use the --mq_txt_folder and --meta_input_file parameters at the same time.")
        
        meta_input_df = pd.read_csv(args.meta_input_file, sep='\t')
        meta_input_df.columns = ['mq_txt_folder', 'raw_folder', 'tmt_correction_file'][:len(meta_input_df.columns)]
        
        raw_folders = convert_to_path_list(meta_input_df['raw_folder'])
        mq_txt_folders = convert_to_path_list(meta_input_df['mq_txt_folder'])
        if 'tmt_correction_file' in meta_input_df.columns:
            tmt_correction_files = convert_to_path_list(meta_input_df['tmt_correction_file'])
        else:
            tmt_correction_files = [Path(args.tmt_reporter_correction_file)]*len(raw_folders)
    else:
        if not args.raw_folder:
            logging.error("Missing --raw_folder argument")
        if not args.mq_txt_folder:
            logging.error("Missing --mq_txt_folder argument")

        raw_folders = [Path(args.raw_folder)]
        mq_txt_folders = [Path(args.mq_txt_folder)]
        tmt_correction_files = [Path(args.tmt_reporter_correction_file)]
    
    return raw_folders, mq_txt_folders, tmt_correction_files


def convert_to_path_list(s):
    return list(map(Path, s.tolist()))


def parse_stringencies(stringencies):
    if stringencies == '':
        stringencies = [20, 15, 10]
    else:
        try:
            stringencies = [int(i) for i in stringencies.split(',')]
        except ValueError:
            logger.error('This is not a stringency list. Please input the list in the following format: 30,25,20,15,10,5')

    return stringencies


if __name__ == '__main__':
    raise NotImplementedError('Do not run this script.')
