import logging
import multiprocessing
from pathlib import Path
import argparse


logger = logging.getLogger(__name__)


class ArgumentParserWithLogger(argparse.ArgumentParser):
    def error(self, message):
        logger.error(f"Error parsing input arguments: {message}")
        super().error(message)


def parse_args(argv):
    apars = ArgumentParserWithLogger(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    apars.add_argument('--mq_txt_folder', default=None, metavar="DIR", required=True,
                       help='''Path to MaxQuant combined/txt output folder.''')

    apars.add_argument('--raw_folder', default=None, metavar="DIR", required=True,
                       help='''Full path to folder containing .raw or .mzML files.''')

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

    # ------------------------------------------------
    args = apars.parse_args(argv)

    return Path(args.mq_txt_folder), Path(args.raw_folder), parse_stringencies(args.stringencies), Path(
        args.output_folder), args.num_threads, Path(args.tmt_reporter_correction_file), args.tmt_ms_level, args.tmt_requantify


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