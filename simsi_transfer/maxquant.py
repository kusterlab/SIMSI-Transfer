import logging
import warnings
from pathlib import Path
from typing import List, Callable, Any

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def process_and_concat(input_folders: List[Any], reading_function: Callable, **kwargs):
    return pd.concat([reading_function(f, **kwargs) for f in input_folders], axis=0)


def read_msmsscans_txt(mq_txt_folder, tmt_requantify):
    """
    Open msms.txt output file and subselect relevant columns
    :param mq_txt_folder: Processing path containing the 'combined' folder from MQ search
    :return: truncated msmsscans.txt dataframe
    """
    cols = ['Raw file', 'Scan number', 'm/z', 'Mass', 'Retention time', 'Precursor full scan number', 'MS scan number']
    if not tmt_requantify:
        cols += [f'Reporter intensity {i}' for i in range(1, 12)]
        cols += [f'Reporter intensity corrected {i}' for i in range(1, 12)]

    try:
        msmsscans = pd.read_csv(mq_txt_folder / Path('msmsScans.txt'), sep='\t', usecols=cols)
    except ValueError:
        cols.remove('Reporter intensity 11')
        cols.remove('Reporter intensity corrected 11')
        msmsscans = pd.read_csv(mq_txt_folder / Path('msmsScans.txt'), sep='\t', usecols=cols)
    
    msmsscans = msmsscans.rename(
            columns={'Scan number': 'scanID'})
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
