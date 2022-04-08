import pandas as pd


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


if __name__ == '__main__':
    raise NotImplementedError('Do not run this script.')


def merge_summary_with_evidence(summary: pd.DataFrame, evidence: pd.DataFrame):
    '''
    take subset of evidence.txt columns that does not overlap with the columns in the summary dataframe
    '''
    merge_columns = ['evidence_ID', 'Modified sequence', 'Raw file', 'Charge', 'Leading proteins', 'Type', 'Calibrated retention time',
                     'Calibrated retention time start', 'Calibrated retention time finish', 'Retention time calibration', 'Intensity']
    evidence = evidence[merge_columns]
    summary = pd.merge(left=summary, right=evidence, left_on=['Modified sequence', 'Raw file', 'Charge'],
                       right_on=['Modified sequence', 'Raw file', 'Charge'], how='left')
    return summary