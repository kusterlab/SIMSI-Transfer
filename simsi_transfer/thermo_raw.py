import os
from sys import platform
from typing import Optional
from pathlib import Path
import logging

from .utils import subprocess_with_logger as subprocess

logger = logging.getLogger(__name__)


def convert_raw_mzml(input_path: Path, output_path: Optional[Path] = None, gzip = False, ms_level: str = "2-"):
    """
    Converts a ThermoRaw file to mzML. Adapted from prosit_io/raw/thermo_raw.py

    :param input_path: File path of the Thermo Rawfile
    :param output_path: File path of the mzML path
    :param ms_level: MS levels to keep, "2-" means MS2 and above (e.g. MS3)
    """
    if output_path is None:
        output_path = input_path.stem + ".mzML"
    
    if output_path.is_file():
        logger.info(f"Found converted file at {output_path}, skipping conversion")
        return output_path
    
    if gzip:
        gzip = "-g"
    else:
        gzip = ""

    if "linux" in platform:
        mono = "mono"
    elif "win" in platform:
        mono = ""
    
    exec_path = Path(__file__).parent.absolute() # get path of parent directory of current file
    exec_command = f"{mono} {exec_path}/utils/ThermoRawFileParser/ThermoRawFileParser.exe {gzip} --msLevel \"{ms_level}\" -i \"{input_path}\" -b \"{output_path}.tmp\""
    logger.info(f"Converting thermo rawfile to mzml with the command: '{exec_command}'")
    subprocess.run(exec_command)
    
    # only rename the file now, so that we don't have a partially converted file if something fails
    os.rename(f"{output_path}.tmp", output_path)
    
    return output_path


def convert_raw_mzml_batch(raw_folder: Path, output_folder: Optional[Path] = None, num_threads: int = 1, gzip = False, ms_level: str = "2-"):
    raw_files = get_raw_files(raw_folder)
    
    if not output_folder.is_dir():
        output_folder.mkdir(parents=True)
    
    if num_threads > 1:
        import multiprocessing
        from .utils.multiprocessing_pool import JobPool
        processingPool = JobPool(processes=num_threads)
    
    mzml_files = []
    for raw_file in raw_files:
        if raw_file.suffix == ".raw":
            input_path = raw_file
            output_path = output_folder / raw_file.with_suffix('.mzML').name
            
            if num_threads > 1:
                processingPool.applyAsync(convert_raw_mzml, (input_path, output_path))
            else:
                raw_file = convert_raw_mzml(input_path, output_path)
                mzml_files.append(raw_file)
    
    if num_threads > 1:
        mzml_files = processingPool.checkPool(printProgressEvery = 1)
        
    return mzml_files
    

def get_raw_files(raw_folder: Path):
    """
    Obtains raw files by scanning through the raw_path directory.
    """
    raw_files = []
    if not raw_folder.is_dir():
        raise ValueError(f'Failed converting raw files, {raw_folder} is not a directory.')
    
    extension = ".raw"
    raw_files = [f for f in raw_folder.iterdir() if f.suffix.lower() == extension]
    
    if len(raw_files) == 0:
        extension = ".mzml"
        raw_files = [f for f in raw_folder.iterdir() if f.suffix.lower() == extension]
    
    if len(raw_files) == 0:
        raise ValueError(f'Failed converting raw files, {raw_folder} did not contain any .mzML or .raw files.')
        
    logger.info(f"Found {len(raw_files)} raw files in the search directory")
    return raw_files
            


if __name__ == "__main__":
    from sys import argv
    if len(argv) == 2:
        converter = convert_raw_mzml(argv[1], ms_level = "2-")
    else:
        logger.error("Please specify a rawfile")
