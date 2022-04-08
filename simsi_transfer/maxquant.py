import logging
from pathlib import Path

import pandas as pd


logger = logging.getLogger(__name__)


def read_msmsscans_txt(mainpath, tmt_requantify):
    """
    Open msms.txt output file and subselect relevant columns
    :param mainpath: Processing path containing the 'combined' folder from MQ search
    :return: truncated msmsscans.txt dataframe
    """
    cols = ['Raw file', 'Scan number', 'm/z', 'Mass', 'Retention time', 'Precursor full scan number', 'MS scan number']
    if not tmt_requantify:
        cols += [f'Reporter intensity {i}' for i in range(1,12)]
        cols += [f'Reporter intensity corrected {i}' for i in range(1,12)]

    try:
        msmsscans = pd.read_csv(mainpath / Path('msmsScans.txt'), sep='\t', usecols=cols).rename(
            columns={'Scan number': 'scanID'})
        tmt = 11
    except ValueError:
        cols.remove('Reporter intensity 11')
        cols.remove('Reporter intensity corrected 11')
        msmsscans = pd.read_csv(mainpath / Path('msmsScans.txt'), sep='\t', usecols=cols).rename(
            columns={'Scan number': 'scanID'})
        tmt = 10
    return msmsscans, tmt


def read_msms_txt(mainpath):
    """
    Open msms.txt output file and subselect relevant columns
    :param mainpath: Processing path containing the 'combined' folder from MQ search
    :return: truncated msms.txt dataframe
    """
    msmstxt = pd.read_csv(mainpath / Path('msms.txt'), sep='\t',
                          usecols=['Raw file', 'Scan number', 'Sequence', 'Modified sequence', 'Length',
                                   'Modifications', 'Missed cleavages', 'Proteins', 'Gene Names', 'Protein Names',
                                   'Charge', 'Mass error [ppm]', 'PIF', 'Precursor Intensity', 'PEP', 'Score',
                                   'Delta score', 'Reverse']).rename(
        columns={'Scan number': 'scanID'})
    return msmstxt


def read_evidence_txt(mainpath):
    """
    Open msms.txt output file and subselect relevant columns
    :param mainpath: Processing path containing the 'combined' folder from MQ search
    :return: truncated evidence.txt dataframe
    """
    usecols = ['Sequence', 'Modified sequence', 'Leading proteins', 'Raw file', 'Experiment', 'Fraction',
               'Charge', 'Calibrated retention time', 'Retention time', 'Retention length',
               'Calibrated retention time start', 'Calibrated retention time finish',
               'Retention time calibration', 'Type', 'Intensity', 'Reverse']
    evidence = pd.read_csv(mainpath / Path('evidence.txt'), sep='\t', usecols=lambda x: x in usecols)
    return evidence


def read_allpeptides_txt(mainpath):
    """
    Open msms.txt output file and subselect relevant columns
    :param mainpath: Processing path containing the 'combined' folder from MQ search
    :return: truncated evidence.txt dataframe
    """
    cols = ['Raw file', 'Type', 'Charge', 'm/z', 'Retention time', 'Retention length', 'Min scan number', 'Max scan number', 'Intensity']
    evidence = pd.read_csv(mainpath / Path('allPeptides.txt'), sep='\t',
                           usecols=cols)
    return evidence


def get_rawfile_metadata(evidence_txt):
    if 'Fraction' in evidence_txt.columns:
        return evidence_txt[['Raw file', 'Experiment', 'Fraction']].drop_duplicates()
    else:
        meta_df = evidence_txt[['Raw file', 'Experiment']].drop_duplicates()
        meta_df['Fraction'] = 1
        return meta_df


if __name__ == '__main__':
    raise NotImplementedError('Do not run this script.')
