import logging
import re
import warnings
from pathlib import Path
from typing import List, Callable, Any

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def process_and_concat(input_folders: List[Any], reading_function: Callable, **kwargs):
    return pd.concat([reading_function(f, **kwargs) for f in input_folders], axis=0)


def get_plex(input_folders):
    columns = pd.read_csv(input_folders[0]/Path('msms.txt'), nrows=0, sep='\t').columns.tolist()
    substring = re.compile(r'^Reporter intensity (\d{1,2})$')
    reporters = [i for i in columns if re.match(substring, i)]
    plex_number = max([int(re.search(substring, i).group(1)) for i in reporters])
    logger.info(f'Automatic plex detection: {plex_number}plex detected.')
    return plex_number


def read_msmsscans_txt(mq_txt_folder, tmt_requantify, plex):
    """
    Open msms.txt output file and subselect relevant columns
    :param mq_txt_folder: Processing path containing the 'combined' folder from MQ search
    :return: truncated msmsscans.txt dataframe
    """
    cols = ['Raw file', 'Scan number', 'm/z', 'Mass', 'Retention time', 'Precursor full scan number', 'MS scan number']
    if not tmt_requantify:
        cols += [f'Reporter intensity {i}' for i in range(1, plex + 1)]
        cols += [f'Reporter intensity corrected {i}' for i in range(1, plex + 1)]
    msmsscans = pd.read_csv(mq_txt_folder / Path('msmsScans.txt'), sep='\t', usecols=cols)
    
    msmsscans = msmsscans.rename(columns={'Scan number': 'scanID'})
    return msmsscans


def read_msms_txt(mq_txt_folder):
    """
    Open msms.txt output file and subselect relevant columns
    :param mq_txt_folder: Processing path containing the 'combined' folder from MQ search
    :return: truncated msms.txt dataframe
    """
    columns = ['Raw file', 'Scan number', 'Sequence', 'Modified sequence', 'Phospho (STY) Probabilities', 'Length',
               'Modifications', 'Missed cleavages', 'Proteins', 'Gene Names', 'Protein Names',
               'Charge', 'Mass error [ppm]', 'PIF', 'Precursor Intensity', 'PEP', 'Score',
               'Delta score', 'Reverse']
    msmstxt = pd.read_csv(mq_txt_folder / Path('msms.txt'), sep='\t',
                          usecols=lambda x: x in columns)
    
    for col in columns:
        if col not in msmstxt.columns:
            logger.warning(f'Missing column in msms.txt, filled with numpy NaN: {col}')
            msmstxt[col] = np.NaN
    
    msmstxt = msmstxt.rename(columns={'Scan number': 'scanID'})
    return msmstxt


def read_evidence_txt(mq_txt_folder):
    """
    Open msms.txt output file and subselect relevant columns
    :param mq_txt_folder: Processing path containing the 'combined' folder from MQ search
    :return: truncated evidence.txt dataframe
    """
    usecols = ['Sequence', 'Modified sequence', 'Leading proteins', 'Raw file', 'Experiment', 'Fraction',
               'Charge', 'Calibrated retention time', 'Retention time', 'Retention length',
               'Calibrated retention time start', 'Calibrated retention time finish',
               'Retention time calibration', 'Type', 'Intensity', 'Reverse']
    evidence = pd.read_csv(mq_txt_folder / Path('evidence.txt'), sep='\t', usecols=lambda x: x in usecols)
    return evidence


def read_allpeptides_txt(mq_txt_folder):
    """
    Open msms.txt output file and subselect relevant columns
    :param mq_txt_folder: Processing path containing the 'combined' folder from MQ search
    :return: truncated evidence.txt dataframe
    """
    cols = ['Raw file', 'Type', 'Charge', 'm/z', 'Retention time', 'Retention length', 'Min scan number',
            'Max scan number', 'Intensity']
    evidence = pd.read_csv(mq_txt_folder / Path('allPeptides.txt'), sep='\t',
                           usecols=cols)
    return evidence


def get_rawfile_metadata(evidence_txt):
    metadata_columns = ['Raw file', 'Experiment', 'Fraction']
    metadata_columns_available = [x for x in metadata_columns if x in evidence_txt.columns]
    metadata_columns_unavailable = [x for x in metadata_columns if x not in evidence_txt.columns]

    meta_df = evidence_txt[metadata_columns_available].drop_duplicates()
    meta_df[metadata_columns_unavailable] = 1
    return meta_df


if __name__ == '__main__':
    raise NotImplementedError('Do not run this script.')
