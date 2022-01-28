import subprocess
from typing import List
from pathlib import Path


def cluster_mzml_files(mzml_files: List[Path], pvals: List[float], maracluster_folder: Path, num_threads: int = 1, cleanup: bool = True):
    """
    Runs maracluster on a list of mzML files
    """
    if not maracluster_folder.is_dir():
        maracluster_folder.mkdir(parents=True)
        
    batch_file = create_batch_file(maracluster_folder, mzml_files)
    
    pvals_string = ",".join(map(str, pvals))
    cluster_command = f'OMP_NUM_THREADS={num_threads} maracluster batch -b {batch_file} -c {pvals_string} -f {maracluster_folder} 2>&1'
    process = subprocess.run(
        cluster_command,
        shell=True,
        check=True)


def create_batch_file(maracluster_folder: Path, mzml_files: List[Path]):
    batch_file = maracluster_folder / Path('file_list.txt')
    with open(batch_file, 'w') as w:
        for mzml_file in mzml_files:
            w.write(str(mzml_file) + "\n")
    return batch_file
