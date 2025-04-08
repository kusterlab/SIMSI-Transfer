from sys import platform
import subprocess
import logging
from typing import List
from pathlib import Path

import pandas as pd

from .utils import subprocess_with_logger as subprocess

logger = logging.getLogger(__name__)


def has_previous_run(mainpath: Path, mzml_files: List[Path], pvals: List[float]):
    if not (mainpath / Path(f'file_list.txt')).is_file():
        return False
    
    with open(mainpath / Path(f'file_list.txt'), 'r') as file:
        file_list = set(Path(line.rstrip()) for line in file)
        mzml_file_list = set(mzml_files)
        if mzml_file_list != file_list:
            logger.warning("Found previous MaRaCluster run with a different file list")
            logger.warning(f"- Missing in file_list.txt: {mzml_file_list.difference(file_list)}")
            logger.warning(f"- Missing in mzml_files: {file_list.difference(mzml_file_list)}")
            logger.warning("Rerunning MaRaCluster")
            return False

    for pval in pvals:
        if not (mainpath / Path(f'MaRaCluster.clusters_p{pval}.tsv')).is_file():
            return False
    
    return True


def read_cluster_results(mainpath, pval):
    maracluster_df = pd.read_csv(mainpath / Path(f'MaRaCluster.clusters_{pval}.tsv'),
                                 sep='\t', names=['Raw file', 'scanID', 'clusterID'], engine="pyarrow")
    maracluster_df['Raw file'] = maracluster_df['Raw file'].apply(get_file_name)
    return maracluster_df


def get_file_name(raw_file):
    """
    Removes full path and file extension from file name in MaRaCluster column
    :param raw_file: Path string of mzML file
    :return: raw file name without path or extension
    """
    return Path(raw_file).stem


def cluster_mzml_files(mzml_files: List[Path], pvals: List[float], maracluster_folder: Path, dat_folder: Path, num_threads: int = 1, num_threads_per_precursor_bin: int = 0):
    """
    Runs maracluster on a list of mzML files
    """
    if not maracluster_folder.is_dir():
        maracluster_folder.mkdir(parents=True)

    batch_file = create_batch_file(maracluster_folder, mzml_files)
    
    if num_threads_per_precursor_bin > 0:
        from job_pool import JobPool
        
        run_maracluster_index(batch_file, maracluster_folder, dat_folder, num_threads)

        dat_bin_files = read_dat_bin_file_list(maracluster_folder)
        
        processingPool = JobPool(
            processes=max([1, num_threads // num_threads_per_precursor_bin]), timeout=10800, write_progress_to_logger=True
        )  # 10,8000 seconds = 3 hours
        for dat_bin_file in dat_bin_files:
            processingPool.applyAsync(run_maracluster_pvalue_precursor_bin, (maracluster_folder, dat_bin_file, num_threads_per_precursor_bin))
        processingPool.checkPool(printProgressEvery=1)

    run_maracluster_batch(batch_file, pvals, maracluster_folder, dat_folder, num_threads)


def read_dat_bin_file_list(maracluster_folder: Path):
    with open(maracluster_folder / "MaRaCluster.dat_file_list.txt", "r") as f:
        return [Path(line[:-1]).name for line in f.readlines()]


def run_maracluster_index(batch_file: Path, maracluster_folder: Path, dat_folder: Path, num_threads: int = 1):        
    run_maracluster(f'index --batch "{batch_file}" --output-folder "{maracluster_folder}" --dat-folder "{dat_folder}"', num_threads=num_threads)


def run_maracluster_pvalue_precursor_bin(maracluster_folder: Path, dat_bin_file: str, num_threads: int = 1):
    """Computes p-values for a dat file corresponding to a precursor m/z range, e.g. 356.dat

    Args:
        maracluster_folder (Path): _description_
        dat_bin_file (str): _description_
        num_threads (int, optional): _description_. Defaults to 1.
    """    
    run_maracluster(f'pvalue --output-folder "{maracluster_folder}" --prefix "{dat_bin_file}" --specIn "{maracluster_folder}/{dat_bin_file}" --peakCountsFN "{maracluster_folder}/MaRaCluster.peak_counts.dat" --clusteringTree "{maracluster_folder}/{dat_bin_file}.pvalue_tree.tsv"', num_threads=num_threads)


def run_maracluster_batch(batch_file: Path, pvals: List[float], maracluster_folder: Path, dat_folder: Path, num_threads: int = 1):
    """
    Runs maracluster on a list of mzML files
    """       
    pvals_string = ",".join(map(lambda x : str(-1*x), pvals))
    run_maracluster(f'batch --batch "{batch_file}" --clusterThresholds {pvals_string} --output-folder "{maracluster_folder}" --dat-folder "{dat_folder}"', num_threads=num_threads)


def run_maracluster(maracluster_cmd: str, num_threads: int = 1):
    """
    Runs maracluster on a list of mzML files
    """
    num_threads_string = f"OMP_NUM_THREADS={num_threads}"
    if "win" in platform:
        num_threads_string = f"set OMP_NUM_THREADS={num_threads} &&"
        
    exec_path = Path(__file__).parent.absolute() # get path of parent directory of current file
    exec_bin = f"{exec_path}/utils/maracluster/linux64/maracluster"
    if "win" in platform:
        exec_bin = f"{exec_path}\\utils\\maracluster\\win64\\maracluster"
        
    cluster_command = f'{num_threads_string} "{exec_bin}" {maracluster_cmd} 2>&1'
    process = subprocess.run(cluster_command)


def create_batch_file(maracluster_folder: Path, mzml_files: List[Path]):
    batch_file = maracluster_folder / Path('file_list.txt')
    with open(batch_file, 'w') as w:
        for mzml_file in mzml_files:
            w.write(str(mzml_file) + "\n")
    return batch_file
