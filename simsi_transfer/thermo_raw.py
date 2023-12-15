import os
from sys import platform
from typing import Optional, List
from pathlib import Path
import logging

from .utils import subprocess_with_logger as subprocess

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + __file__)


def convert_raw_mzml(input_path: Path, output_path: Optional[Path] = None, gzip: bool = False, ms_level: str = "2-") -> Path:
    """
    Converts a ThermoRaw file to mzML. Adapted from prosit_io/raw/thermo_raw.py

    :param input_path: File path of the Thermo Rawfile
    :param output_path: File path of the mzML path
    :param ms_level: MS levels to keep, "2-" means MS2 and above (e.g. MS2 and MS3)
    """
    if output_path is None:
        output_path = input_path.stem + ".mzML"
    
    if output_path.is_file():
        logger.debug(f"Found converted file at {output_path}, skipping conversion")
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
    logger.debug(f"Converting thermo rawfile to mzml with the command: '{exec_command}'")
    subprocess.run(exec_command)
    
    # only rename the file now, so that we don't have a partially converted file if something fails
    os.rename(f"{output_path}.tmp", output_path)
    
    return output_path


def convert_raw_mzml_batch(raw_files: List[Path], output_folder: Optional[Path] = None, num_threads: int = 1, gzip: bool = False, ms_level: str = "2-") -> List[Path]:
    if not output_folder.is_dir():
        output_folder.mkdir(parents=True)
    
    if num_threads > 1:
        from job_pool import JobPool
        processingPool = JobPool(processes=num_threads, write_progress_to_logger=True)
    
    mzml_files = []
    for raw_file in raw_files:
        if raw_file.suffix.lower() == ".raw":
            input_path = raw_file
            output_path = output_folder / raw_file.with_suffix('.mzML').name
            
            if num_threads > 1:
                processingPool.applyAsync(convert_raw_mzml, (input_path, output_path))
            else:
                mzml_file = convert_raw_mzml(input_path, output_path)
                mzml_files.append(mzml_file)
        elif raw_file.suffix.lower() == ".mzml":
            mzml_files.append(raw_file)
    
    if len(mzml_files) == 0 and num_threads > 1:
        mzml_files = processingPool.checkPool()
        
    return mzml_files
    

def get_raw_files(raw_folder: str) -> List[Path]:
    """
    Obtains raw files by scanning through the raw_path directory.
    """
    raw_folder = Path(raw_folder)
    
    raw_files = []
    if not raw_folder.is_dir():
        raise ValueError(f'Failed converting raw files, {raw_folder} is not a directory.')
    
    extension = ".raw"
    raw_files = [f for f in raw_folder.iterdir() if f.suffix.lower() == extension]
    
    if len(raw_files) == 0:
        extension = ".mzml"
        raw_files = [f for f in raw_folder.iterdir() if f.suffix.lower() == extension]
    
    if len(raw_files) == 0:
        raise ValueError(f'Failed getting raw files, {raw_folder} did not contain any .mzML or .raw files.')

    raw_files.sort()
    logger.debug(f"Found {len(raw_files)} raw files in the search directory")
    return raw_files


def check_valid_mzml(mzml_file: Path):
    from pyteomics import mzml
    import lxml

    logger.info(f"Checking if {mzml_file} can be read by pyteomics")
    try:
        with mzml.read(str(mzml_file)) as reader:
            ids = list()
            for i, item in enumerate(reader):
                ids.append((i, item['id']))
    except lxml.etree.XMLSyntaxError as e:
        logger.error("Found invalid symbol, check the error message for the line in the mzML where parsing failed.")
        raise e

    logger.info(f"{mzml_file} was successfully read by pyteomics")


if __name__ == "__main__":
    from sys import argv
    import argparse
    from .command_line_interface import ArgumentParserWithLogger
    from . import __version__, __copyright__

    desc = f'SIMSI-raw-convert version {__version__}\n{__copyright__}' 
    apars = ArgumentParserWithLogger(
        description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    apars.add_argument('--raw_file', default=None, metavar="R", required=True,
                       help='''
                       Path to raw file you want to convert.
                       ''')

    apars.add_argument('--mzml_output_file', default=None, metavar="M",
                       help='''Output path to mzML file.''')
    
    apars.add_argument('--validate', default=False, action='store_true',
                       help='''Check if mzML file can be read by pyteomics.''')
    
    args = apars.parse_args(argv[1:])

    logger.info(f'{desc}')
    logger.info(f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}')

    converter = convert_raw_mzml(Path(args.raw_file), Path(args.mzml_output_file), ms_level = "2-")

    if args.validate:
        check_valid_mzml(Path(args.mzml_output_file))