import os
import pandas as pd
from pathlib import Path


def merge_with_msmstxt(tempfile, msms_df):
    """
    Merging function to generate merged dataframe from MaRaCluster clustering results and MaxQuant msms.txt output file,
    requires correct folder structure to work
    :param tempfile: Temporary processing file, e.g. MaRaCluster clustering dataframe
    :param msms_df: MaxQuant msms.txt dataframe
    :return: Merged dataframe
    """
    return pd.merge(left=tempfile, right=msms_df, on=['Raw file', 'scanID'], how='left')


def merge_with_msmsscanstxt(tempfile, msmsscans_df):
    """
    Merging function to generate merged dataframe from MaRaCluster clustering results and MaxQuant msms.txt output file,
    requires correct folder structure to work
    :param tempfile: Temporary processing file, e.g. MaRaCluster clustering dataframe
    :param msmsscans_df: MaxQuant msmsscans.txt dataframe
    :return: Merged dataframe
    """
    return pd.merge(left=tempfile, right=msmsscans_df, on=['Raw file', 'scanID'], how='left')


def merge_with_summarytxt(tempfile, summary_df):
    """
    Merging function to get experiment ID from MaxQuant summary.txt
    :param tempfile:
    :param summary_df:
    :return:
    """
    return pd.merge(left=tempfile, right=summary_df, on='Raw file', how='left')


def assemble_corrected_tmt_table(extracted_folder):
    corrected_tmt = pd.DataFrame(
        columns=['raw_file', 'scanID', 'raw_TMT1', 'raw_TMT2', 'raw_TMT3', 'raw_TMT4', 'raw_TMT5', 'raw_TMT6',
                 'raw_TMT7', 'raw_TMT8', 'raw_TMT9', 'raw_TMT10', 'raw_TMT11', 'raw_TMT12', 'raw_TMT13',
                 'corr_TMT1',
                 'corr_TMT2', 'corr_TMT3', 'corr_TMT4', 'corr_TMT5', 'corr_TMT6', 'corr_TMT7', 'corr_TMT8',
                 'corr_TMT9', 'corr_TMT10', 'corr_TMT11'])

    for file in os.listdir(extracted_folder):
        file = pd.read_csv(Path(extracted_folder / file), sep='\t')
        corrected_tmt = corrected_tmt.append(file)
    corrected_tmt = corrected_tmt.reset_index(drop=True)
    corrected_tmt = corrected_tmt.rename(columns={
        'raw_file': 'Raw file',
        'raw_TMT1': 'Reporter intensity 1',
        'raw_TMT2': 'Reporter intensity 2',
        'raw_TMT3': 'Reporter intensity 3',
        'raw_TMT4': 'Reporter intensity 4',
        'raw_TMT5': 'Reporter intensity 5',
        'raw_TMT6': 'Reporter intensity 6',
        'raw_TMT7': 'Reporter intensity 7',
        'raw_TMT8': 'Reporter intensity 8',
        'raw_TMT9': 'Reporter intensity 9',
        'raw_TMT10': 'Reporter intensity 10',
        'raw_TMT11': 'Reporter intensity 11',
        'corr_TMT1': 'Reporter intensity corrected 1',
        'corr_TMT2': 'Reporter intensity corrected 2',
        'corr_TMT3': 'Reporter intensity corrected 3',
        'corr_TMT4': 'Reporter intensity corrected 4',
        'corr_TMT5': 'Reporter intensity corrected 5',
        'corr_TMT6': 'Reporter intensity corrected 6',
        'corr_TMT7': 'Reporter intensity corrected 7',
        'corr_TMT8': 'Reporter intensity corrected 8',
        'corr_TMT9': 'Reporter intensity corrected 9',
        'corr_TMT10': 'Reporter intensity corrected 10',
        'corr_TMT11': 'Reporter intensity corrected 11'})
    corrected_tmt['Raw file'] = corrected_tmt['Raw file'].str.replace(r'.mzML$', '')
    return corrected_tmt


def merge_with_corrected_tmt(tempfile, corrected_tmt):
    return pd.merge(left=tempfile, right=corrected_tmt, on=['Raw file', 'scanID'], how='left', validate='one_to_one')


if __name__ == '__main__':
    raise NotImplementedError('Do not run this script.')
