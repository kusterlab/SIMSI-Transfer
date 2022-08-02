import time
import os
import sys
import re
from typing import List
from pathlib import Path
import logging

import pandas as pd
import numpy as np
from pyteomics import mzml

logger = logging.getLogger(__name__)

TMT_RAW_COLUMNS = [f'raw_TMT{i}' for i in range(1, 14)]
TMT_CORRECTED_COLUMNS = [f'corr_TMT{i}' for i in range(1, 12)]


def get_correction_factors(correction_factor_path: Path):
    correction = np.array([[100, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 126 C Tag
                           [0.0, 100, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 127 N Tag
                           [0.0, 0.0, 100, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 127 C Tag
                           [0.0, 0.0, 0.0, 100, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 128 N Tag
                           [0.0, 0.0, 0.0, 0.0, 100, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 128 C Tag
                           [0.0, 0.0, 0.0, 0.0, 0.0, 100, 0.0, 0.0, 0.0, 0.0, 0.0],  # 129 N Tag
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100, 0.0, 0.0, 0.0, 0.0],  # 129 C Tag
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100, 0.0, 0.0, 0.0],  # 130 N Tag
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100, 0.0, 0.0],  # 130 C Tag
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100, 0.0],  # 131 N Tag
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100],  # 131 C Tag
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 132 N Overflow
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # 132 C Overflow
                           ])

    # Theoretical TMT Masses in m/z
    tmt_masses = np.array([126.127726, 127.124761, 127.131081, 128.128116, 128.134436, 129.131471,
                           129.137790, 130.134825, 130.141145, 131.138180, 131.144499, 132.141535, 132.147854])


    if correction_factor_path.is_file():
        correction_dataframe = pd.read_csv(correction_factor_path, sep='\t')

        for i in range(11):
            if i not in [0, 1, 2, 3]:
                correction[i - 4, i] = correction_dataframe.iloc[i]['Correction factor -2 [%]']
            if i not in [0, 1]:
                correction[i - 2, i] = correction_dataframe.iloc[i]['Correction factor -1 [%]']
            correction[i + 2, i] = correction_dataframe.iloc[i]['Correction factor +1 [%]']
            if i not in [9, 10]:
                correction[i + 4, i] = correction_dataframe.iloc[i]['Correction factor +2 [%]']

    # Normalize correction factors
    correction_normalized = (correction / correction.sum(axis=0))

    return tmt_masses, correction_normalized


def extract_tmt_reporters(mzml_files: List[Path], output_path: Path, correction_factor_path: Path, num_threads: int = 1,
                          extraction_level: int = 3):
    """
    Takes about 1.5 minute for a 700MB file with 40k MS2 scans
    """
    if not output_path.is_dir():
        output_path.mkdir(parents=True)

    if num_threads > 1:
        from .utils.multiprocessing_pool import JobPool
        processing_pool = JobPool(processes=num_threads)
    for mzml_file in mzml_files:
        if num_threads > 1:
            processing_pool.applyAsync(extract_and_correct_reporters,
                                       (mzml_file, output_path, correction_factor_path, extraction_level))
        else:
            extract_and_correct_reporters(mzml_file, output_path, correction_factor_path, extraction_level)

    if num_threads > 1:
        processing_pool.checkPool(printProgressEvery=1)


def extract_and_correct_reporters(mzml_file, output_path, correction_factor_path, extraction_level):
    tmt_masses, correction_normalized = get_correction_factors(correction_factor_path)

    tolerance = 6 * 1e-3 / 2
    tmt_upper = tmt_masses + tolerance
    tmt_lower = tmt_masses - tolerance

    convert_dict_raw = {k: 'float64' for k in TMT_RAW_COLUMNS}
    convert_dict_corr = {k: 'float64' for k in TMT_CORRECTED_COLUMNS}
    convert_dict_other = {'raw_file': 'str', 'scanID': 'int'}
    convert_dict = {**convert_dict_raw, **convert_dict_corr, **convert_dict_other}
    dfcol = convert_dict.keys()

    output_file = get_extracted_tmt_file_name(output_path, mzml_file)
    if Path(output_file).is_file():
        logger.info(f"Found extracted reporter ions at {output_file}, skipping extraction")
        return

    logger.info('Performing extraction for ' + mzml_file.name)
    fileframe = pd.DataFrame(columns=dfcol)

    seen_scan_ids = set()
    with mzml.read(str(mzml_file)) as reader:
        for i, item in enumerate(reader):
            if item['ms level'] != extraction_level:
                continue

            scanseries = pd.Series(index=dfcol, dtype='float64')

            if extraction_level == 2:
                scanseries['scanID'] = re.search(r'scan=(\d+)', item['id'])[1]
            else:
                # supposed to find parent MS2 spectrum for MS3 by looking into precursorList/precursor/spectrumRef
                scanseries['scanID'] = re.search(r'scan=(\d+)', item['precursorList']['precursor'][0]['spectrumRef'])[1]

            if scanseries['scanID'] in seen_scan_ids:
                logger.warning(
                    f"Found duplicate MS3 spectrum for MS2 spectrum with scan number {scanseries['scanID']} in {mzml_file}, known bug in ThermoRawFileParser...")
                continue
            seen_scan_ids.add(scanseries['scanID'])

            mz = np.array(item['m/z array'])
            intensity = np.array(item['intensity array'])
            for c, (low, upp) in enumerate(zip(tmt_lower, tmt_upper)):
                start_idx = int(np.searchsorted(mz, low))
                end_idx = int(np.searchsorted(mz, upp))
                scanseries[f'raw_TMT{c + 1}'] = intensity[start_idx:end_idx].sum()
            fileframe = pd.concat([fileframe, scanseries.to_frame().T], ignore_index=True)
    fileframe['raw_file'] = mzml_file.name

    fileframe = fileframe.astype(convert_dict)

    # TMT correction
    logger.info('Extraction done, correcting TMT reporters for ' + mzml_file.name)
    fileframe[TMT_CORRECTED_COLUMNS] = pd.DataFrame(fileframe[TMT_RAW_COLUMNS].apply(
        lambda tmt: np.linalg.lstsq(correction_normalized, tmt, rcond=None)[0].round(2), axis=1).tolist(),
                                           columns=TMT_CORRECTED_COLUMNS, index=fileframe[TMT_CORRECTED_COLUMNS].index)
    fileframe[TMT_CORRECTED_COLUMNS] = fileframe[TMT_CORRECTED_COLUMNS].where(fileframe[TMT_CORRECTED_COLUMNS] > 10, 0)
    fileframe.to_csv(output_file, sep='\t', index=False)


def get_extracted_tmt_file_name(output_path: str, mzml_file: Path):
    return f'{output_path}/ext_{mzml_file.name}.txt'


def assemble_corrected_tmt_table(mzml_files, extracted_folder):
    logger.info('Assembling corrected reporter ion tables')
    corrected_tmt = pd.DataFrame(
        columns=['raw_file', 'scanID'] + \
                TMT_RAW_COLUMNS + \
                TMT_CORRECTED_COLUMNS)
    
    corrected_tmts = list()
    for mzml_file in mzml_files:
        tmt_file = get_extracted_tmt_file_name(extracted_folder, mzml_file)
        logger.debug(f'  Appending {tmt_file}')
        df = pd.read_csv(tmt_file, sep='\t')
        corrected_tmts.append(df)
    corrected_tmt = pd.concat(corrected_tmts)
    
    corrected_tmt = corrected_tmt.reset_index(drop=True)
    corrected_tmt = corrected_tmt.rename(columns={
        'raw_file': 'Raw file',
        **{f'raw_TMT{i}': f'Reporter intensity {i}' for i in range(1,12)},
        **{f'corr_TMT{i}': f'Reporter intensity corrected {i}' for i in range(1,12)}
        })
    corrected_tmt['Raw file'] = corrected_tmt['Raw file'].str.replace(r'.mzML$', '', regex=True)
    return corrected_tmt


def merge_with_corrected_tmt(msmsscans_df, corrected_tmt):
    logger.info('Merging corrected reporter ion tables into msmsScans.txt')
    return pd.merge(left=msmsscans_df, right=corrected_tmt, on=['Raw file', 'scanID'], how='left', validate='one_to_one')


if __name__ == "__main__":
    input_files_arg = Path(sys.argv[1])
    extraction_level_arg = int(sys.argv[2])
    output_path_arg = Path(sys.argv[3])

    extract_tmt_reporters([input_files_arg], output_path_arg, extraction_level=extraction_level_arg)

