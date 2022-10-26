from sys import platform
import subprocess
import logging
import re
from typing import List
from pathlib import Path

import pandas as pd

from .utils import subprocess_with_logger as subprocess

logger = logging.getLogger(__name__)


def read_cluster_results(mainpath, pval):
    maracluster_df = pd.read_csv(mainpath / Path(f'MaRaCluster.clusters_{pval}.tsv'),
                                 sep='\t', names=['Raw file', 'scanID', 'clusterID'])
    maracluster_df['Raw file'] = maracluster_df['Raw file'].apply(get_file_name)
    return maracluster_df


def get_file_name(raw_file):
    """
    Removes full path and file extension from file name in MaRaCluster column
    :param raw_file: Path string of mzML file
    :return: raw file name without path or extension
    """
    return Path(raw_file).stem


def cluster_mzml_files(mzml_files: List[Path], pvals: List[float], maracluster_folder: Path, dat_folder: Path, num_threads: int = 1, cleanup: bool = True):
    """
    Runs maracluster on a list of mzML files
    """
    if not maracluster_folder.is_dir():
        maracluster_folder.mkdir(parents=True)
        
    batch_file = create_batch_file(maracluster_folder, mzml_files)
    
    pvals_string = ",".join(map(lambda x : str(-1*x), pvals))
    num_threads_string = f"OMP_NUM_THREADS={num_threads}"
    if "win" in platform:
        num_threads_string = f"set OMP_NUM_THREADS={num_threads} &&"
        
    exec_path = Path(__file__).parent.absolute() # get path of parent directory of current file
    exec_bin = f"{exec_path}/utils/maracluster/linux64/maracluster"
    if "win" in platform:
        exec_bin = f"{exec_path}\\utils\\maracluster\\win64\\maracluster"
        
    cluster_command = f'{num_threads_string} "{exec_bin}" batch --batch "{batch_file}" --clusterThresholds {pvals_string} --output-folder "{maracluster_folder}" --dat-folder "{dat_folder}" 2>&1'
    process = subprocess.run(cluster_command)


def create_batch_file(maracluster_folder: Path, mzml_files: List[Path]):
    batch_file = maracluster_folder / Path('file_list.txt')
    with open(batch_file, 'w') as w:
        for mzml_file in mzml_files:
            w.write(str(mzml_file) + "\n")
    return batch_file
