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

pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


def get_correction_factors(mode: str, correction_factor_path: Path):
    if mode == 'topas':
        # TMT 11-plex Correction factors for LOT:
        # Reporter percentage   26   27N  27C  28N  28C  29N  29C  30N  30C  31N
        correction = np.array([[100, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 126 C Tag
                               [0.0, 100, 0.0, 0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 127 N Tag
                               [6.9, 0.0, 100, 0.0, 1.4, 0.0, 0.0, 0.0, 0.0, 0.0],  # 127 C Tag
                               [0.0, 7.3, 0.0, 100, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0],  # 128 N Tag
                               [0.1, 0.0, 5.9, 0.0, 100, 0.0, 2.3, 0.0, 0.0, 0.0],  # 128 C Tag
                               [0.0, 0.2, 0.0, 5.7, 0.0, 100, 0.0, 1.4, 0.0, 0.0],  # 129 N Tag
                               [0.0, 0.0, 0.0, 0.0, 4.9, 0.0, 100, 0.0, 1.7, 0.0],  # 129 C Tag
                               [0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 100, 0.0, 2.1],  # 130 N Tag
                               [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.3, 0.0, 100, 0.0],  # 130 C Tag
                               [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.3, 0.0, 100],  # 131 N Tag
                               [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.8, 0.0],  # 131 C Overflow
                               [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.3, 0.0, 1.7]  # 132 N Overflow
                               ])
        # Theoretical TMT Masses in m/z
        TMT = np.array([126.127726, 127.124761, 127.131081, 128.128116, 128.134436, 129.131471,
                        129.137790, 130.134825, 130.141145, 131.138180, 131.144499, 132.141535])

        tmt_raw_col = ['raw_TMT1', 'raw_TMT2', 'raw_TMT3', 'raw_TMT4', 'raw_TMT5', 'raw_TMT6', 'raw_TMT7', 'raw_TMT8',
                       'raw_TMT9', 'raw_TMT10', 'raw_TMT11', 'raw_TMT12']
        tmt_corr_col = ['corr_TMT1', 'corr_TMT2', 'corr_TMT3', 'corr_TMT4', 'corr_TMT5', 'corr_TMT6', 'corr_TMT7',
                        'corr_TMT8', 'corr_TMT9', 'corr_TMT10']

    elif mode == 'custom':
        # # custom TMT correction factors: 8plex
        # # Reporter percentage   26   27N  27C  28C  29N  29C  30C  31N
        # correction = np.array([[100, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 126 C Tag
        #                        [0.0, 100, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 127 N Tag
        #                        [6.3, 0.0, 100, 1.4, 0.0, 0.0, 0.0, 0.0],  # 127 C Tag
        #                        [0.0, 0.0, 5.7, 100, 0.0, 1.5, 0.0, 0.0],  # 128 C Tag
        #                        [0.0, 0.0, 0.0, 0.0, 100, 0.0, 0.0, 0.0],  # 129 N Tag
        #                        [0.0, 0.0, 0.0, 4.5, 0.0, 100, 3.3, 0.0],  # 129 C Tag
        #                        [0.0, 0.0, 0.0, 0.0, 0.0, 3.7, 100, 0.0],  # 130 C Tag
        #                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100],  # 131 N Tag
        #                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 131 C Overflow
        #                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # 132 N Overflow
        #                        ])
        # Theoretical TMT Masses in m/z
        #
        # TMT = np.array([126.127726, 127.124761, 127.131081, 128.134436, 129.131471,
        #                 129.137790, 130.141145, 131.138180, 131.144499, 132.141535])
        #
        # tmt_raw_col = ['raw_TMT1', 'raw_TMT2', 'raw_TMT3', 'raw_TMT4', 'raw_TMT5', 'raw_TMT6', 'raw_TMT7', 'raw_TMT8',
        #                'raw_TMT9', 'raw_TMT10']
        # tmt_corr_col = ['corr_TMT1', 'corr_TMT2', 'corr_TMT3', 'corr_TMT4', 'corr_TMT5', 'corr_TMT6', 'corr_TMT7',
        #                 'corr_TMT8']

        # # custom TMT correction factors: 10plex
        # # Reporter percentage   26   27N  27C  28N  28C  29N  29C  30N  30C  31N
        # correction = np.array([[100, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 126 C Tag
        #                        [0.0, 100, 0.0, 0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 127 N Tag
        #                        [7.2, 0.0, 100, 0.0, 1.4, 0.0, 0.0, 0.0, 0.0, 0.0],  # 127 C Tag
        #                        [0.0, 7.3, 0.0, 100, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0],  # 128 N Tag
        #                        [0.2, 0.0, 6.3, 0.0, 100, 0.0, 2.3, 0.0, 0.0, 0.0],  # 128 C Tag
        #                        [0.0, 0.2, 0.0, 5.7, 0.0, 100, 0.0, 2.7, 0.0, 0.0],  # 129 N Tag
        #                        [0.0, 0.0, 0.0, 0.0, 5.1, 0.0, 100, 0.0, 2.9, 0.0],  # 129 C Tag
        #                        [0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 100, 0.0, 3.4],  # 130 N Tag
        #                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.3, 0.0, 100, 0.0],  # 130 C Tag
        #                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.9, 0.0, 100],  # 131 N Tag
        #                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.3, 0.0],  # 131 C Overflow
        #                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.3, 0.0, 3.3]  # 132 N Overflow
        #                        ])

        # custom TMT correction factors: 10plex Krug et al
        # Reporter percentage   26   27N  27C  28N  28C  29N  29C  30N  30C  31N
        # correction = np.array([[100, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 126 C Tag
        #                        [0.0, 100, 0.0, 0.4, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0],  # 127 N Tag
        #                        [4.8, 0.0, 100, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0],  # 127 C Tag
        #                        [0.0, 4.7, 0.0, 100, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0],  # 128 N Tag
        #                        [0.0, 0.0, 4.4, 0.0, 100, 0.0, 1.3, 0.0, 0.0, 0.0],  # 128 C Tag
        #                        [0.0, 0.0, 0.0, 3.4, 0.0, 100, 0.0, 1.3, 0.0, 0.3],  # 129 N Tag
        #                        [0.0, 0.0, 0.0, 0.0, 3.8, 0.0, 100, 0.0, 1.8, 0.0],  # 129 C Tag
        #                        [0.0, 0.0, 0.0, 0.0, 0.0, 3.1, 0.0, 100, 0.0, 1.7],  # 130 N Tag
        #                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.9, 0.0, 100, 0.0],  # 130 C Tag
        #                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.8, 0.0, 100],  # 131 N Tag
        #                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.1, 0.0],  # 131 C Overflow
        #                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.7, 0.0, 1.6]  # 132 N Overflow
        #                        ])

        correction_dataframe = pd.read_csv(correction_factor_path, sep='\t')
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
                               [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]   # 132 C Overflow
                               ])

        for i in range(11):
            if i not in [0, 1, 2, 3]:
                correction[i - 4, i] = correction_dataframe.iloc[i]['Correction factor -2 [%]']
            if i not in [0, 1]:
                correction[i - 2, i] = correction_dataframe.iloc[i]['Correction factor -1 [%]']
            correction[i + 2, i] = correction_dataframe.iloc[i]['Correction factor +1 [%]']
            if i not in [9, 10]:
                correction[i + 4, i] = correction_dataframe.iloc[i]['Correction factor +2 [%]']

        # Theoretical TMT Masses in m/z
        TMT = np.array([126.127726, 127.124761, 127.131081, 128.128116, 128.134436, 129.131471,
                        129.137790, 130.134825, 130.141145, 131.138180, 131.144499, 132.141535, 132.147854])

        tmt_raw_col = ['raw_TMT1', 'raw_TMT2', 'raw_TMT3', 'raw_TMT4', 'raw_TMT5', 'raw_TMT6', 'raw_TMT7', 'raw_TMT8',
                       'raw_TMT9', 'raw_TMT10', 'raw_TMT11', 'raw_TMT12', 'raw_TMT13']
        tmt_corr_col = ['corr_TMT1', 'corr_TMT2', 'corr_TMT3', 'corr_TMT4', 'corr_TMT5', 'corr_TMT6', 'corr_TMT7',
                        'corr_TMT8', 'corr_TMT9', 'corr_TMT10', 'corr_TMT11']

    else:
        raise ValueError("Please select 'topas' or 'custom' mode.")

    # Normalize correction factors
    correction_normalized = (correction / correction.sum(axis=0))
    
    return TMT, tmt_raw_col, tmt_corr_col, correction_normalized


def extract_tmt_reporters(mzml_files: List[Path], output_path: Path, correction_factor_path:Path, num_threads: int = 1, mode: str = "custom", extraction_level: int = 3):
    """
    Takes about 1.5 minute for a 700MB file with 40k MS2 scans
    """
    if not output_path.is_dir():
        output_path.mkdir(parents=True)
    
    TMT, tmt_raw_col, tmt_corr_col, correction_normalized = get_correction_factors(mode, correction_factor_path)
    
    tolerance = 6 * 1e-3 / 2
    tmt_upper = TMT + tolerance
    tmt_lower = TMT - tolerance

    dfcol = ['raw_file', 'scanID'] + tmt_raw_col + tmt_corr_col
    
    for mzml_file in mzml_files:
        output_file = f'{output_path}/ext_{mzml_file.name}.txt'
        if Path(output_file).is_file():
            logger.info(f"Found extracted reporter ions at {output_file}, skipping extraction")
            continue
        
        logger.info('Performing extraction for ' + mzml_file.name)
        fileframe = pd.DataFrame(columns=dfcol)
        with mzml.read(str(mzml_file)) as reader:
            for i, item in enumerate(reader):
                if i % 1000 == 0:
                    logger.info(f"Processing spectrum {i}")

                if item['ms level'] != extraction_level:
                    continue

                scanseries = pd.Series(index=dfcol, dtype='float64')

                if extraction_level == 2:
                    scanseries['scanID'] = re.search(r'scan=(\d+)', item['id'])[1]
                else:
                    # supposed to find parent MS2 spectrum for MS3 by looking into precursorList/precursor/spectrumRef
                    scanseries['scanID'] = re.search(r'scan=(\d+)', item['precursorList']['precursor'][0]['spectrumRef'])[1]

                mz = np.array(item['m/z array'])
                intensity = np.array(item['intensity array'])
                for c, (low, upp) in enumerate(zip(tmt_lower, tmt_upper)):
                    start_idx = int(np.searchsorted(mz, low))
                    end_idx = int(np.searchsorted(mz, upp))
                    scanseries[f'raw_TMT{c + 1}'] = intensity[start_idx:end_idx].sum()
                fileframe = fileframe.append(scanseries, ignore_index=True)
        fileframe['raw_file'] = mzml_file.name

        # TMT correction
        logger.info('Extraction done, correcting TMT reporters for ' + mzml_file.name)
        fileframe[tmt_corr_col] = pd.DataFrame(fileframe[tmt_raw_col].apply(
            lambda tmt: np.linalg.lstsq(correction_normalized, tmt, rcond=None)[0].round(2), axis=1).tolist(),
                                               columns=tmt_corr_col, index=fileframe[tmt_corr_col].index)
        fileframe[tmt_corr_col] = fileframe[tmt_corr_col].where(fileframe[tmt_corr_col] > 10, 0)
        fileframe['scanID'] = fileframe['scanID'].astype(int)
        fileframe.to_csv(
            output_file, sep='\t', index=False)


if __name__ == "__main__":
    input_files = Path(sys.argv[1])
    extraction_level = int(sys.argv[2])
    output_path = Path(sys.argv[3])
    
    extract_tmt_reporters([input_files], output_path, extraction_level=extraction_level)
    
