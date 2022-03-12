import re
import math
import os
from sys import platform
import subprocess
from pathlib import Path
import logging
import datetime
from typing import List

import pandas as pd
import numpy as np

from .merging_functions import merge_with_msmsscanstxt, merge_with_msmstxt, merge_with_summarytxt


logger = logging.getLogger(__name__)


def purge_mrc_files(raw_file, mode='mzML'):
    """
    Removes full path and file extension from file name in MaRaCluster column, default mzML files (case insensitive)
    :param raw_file: Path string of mzML file
    :param mode: file extension string, defaults to mzML
    :return: raw file name without path or extension
    """
    rx = r'.+/'
    if "win" in platform:
        rx = r'.+\\'
    
    while True:
        if re.search(rx, raw_file):
            raw_file = re.sub(rx, '', raw_file)
            continue
        break
    raw_file = re.sub(fr'\.{mode}', '', raw_file, flags=re.IGNORECASE)
    return raw_file


def purge_modified_sequence(pepseq):
    """
    Completely removes all modification strings in parentheses, leaving only raw sequences
    :param pepseq: Modified peptide sequence
    :return: Raw peptide sequence
    """
    if type(pepseq) == str:
        rx = re.compile(r'\([^()]+\)')  # remove any annotation in brackets
        while True:
            if re.search(rx, pepseq):
                pepseq = re.sub(rx, '', pepseq)
                continue
            break
        rx = re.compile(r'_')
        if re.search(rx, pepseq):
            pepseq = re.sub(rx, '', pepseq)
        rx = re.compile(r'p([STY])')
        if re.search(rx, pepseq):
            pepseq = re.sub(rx, r'\g<1>', pepseq)
        return pepseq
    elif math.isnan(pepseq):
        return np.nan
    else:
        raise ValueError('String required for modified sequences')


def clean_modified_sequence(pepseq):
    """
    Normalizing modified sequences into p[STY] annotation
    :param pepseq: Modified peptide sequence with X(<Modification>(<Residue>)) annotation and leading or trailing
    underscores
    :return: Modified peptide sequence with p[STY] annotation and no other modifications
    :raises: ValueError if other datatype than string or np.nan is used as input
    """
    if type(pepseq) == str:
        rx = re.compile(r'_')
        if re.search(rx, pepseq):
            pepseq = re.sub(rx, '', pepseq)
        rx = re.compile(r'\(Acetyl \(Protein N-term\)\)')
        if re.search(rx, pepseq):
            pepseq = re.sub(rx, '', pepseq)
        rx = re.compile(r'\(Oxidation \(M\)\)')
        if re.search(rx, pepseq):
            pepseq = re.sub(rx, '', pepseq)
        rx = re.compile(r'([STY])\(Phospho \(STY\)\)')
        if re.search(rx, pepseq):
            pepseq = re.sub(rx, r'p\g<1>', pepseq)
        # TODO: Add conversion from older MQ version, THEPEPT(ph)IDE, +ac +ox
        return pepseq
    elif pepseq != pepseq:
        return np.nan
    else:
        raise ValueError('String required for modified sequences')


def remove_nan_values(input_list):
    """
    Eliminates NaNs from sets or lists, returns list of all other values
    :param input_list: List (or set) of elements
    :return: List with removed np.NaN values
    """
    return [v for v in input_list if pd.notnull(v)]


def get_single_item(input_list):
    """
    Extracts single item from list or set of length = 1 and returns it
    :param input_list: List (or set) of len = 1
    :return: Element of that list
    :raises: IndexError, if initial list has more than one element
    """
    if len(input_list) == 1:
        return next(iter(input_list))
    else:
        raise IndexError('Multiple elements found in input, should only contain one element')


def transfer(sumdf, rawseq='Sequence', modseq='Modified sequence', mask=False):
    """
    Main function for transfers by clustering. Transfers identifications for merged dataframe and adds a column for
    identification type. Transferred columns are Sequence, Modified sequence, Proteins, Gene names, Protein Names,
    Charge, m/z, and Mass.
    :param sumdf: Summary dataframe, merged from cleaned msms.txt and MaRaCluster clusters.tsv file
    :param rawseq: Column name of the raw sequence column
    :param modseq: Column name of the modified sequence column
    :param mask: if false, uses "identification" column from summary frame. If set to a value, transfer uses the
    'identification_{mask}' column, needed for the masking analysis
    :return: DataFrame with transferred identifications resembling MaxQuant msmsScans.txt
    """
    if mask:
        ident = f'identification_{mask}'
    else:
        ident = 'identification'
    grpdf = sumdf[sumdf[modseq] == sumdf[modseq]].groupby('clusterID')[modseq].unique().reset_index(
        name='elements')
    grpdf['setcount'] = grpdf['elements'].str.len()  # set length of sequences without NaN
    grpdf = grpdf[grpdf['setcount'] == 1]  # Only sets with one sequence are allowed
    grpdf[ident] = 't'
    grpdf['representative_transfer_sequence'] = grpdf['elements'].apply(get_single_item)
    grpdf.drop(columns=['elements', 'setcount'], axis=1, inplace=True)
    mg1 = sumdf.merge(grpdf, on=['clusterID'], how='left')
    mg1.loc[mg1[modseq] == mg1[modseq], ident] = 'd'
    mg1.loc[mg1[modseq] != mg1[modseq], modseq] = mg1['representative_transfer_sequence']
    mg1.drop(columns=['representative_transfer_sequence'], axis=1, inplace=True)
    
    replacement_dict = {rawseq: 'd_Sequence', 
                        'Modifications': 'd_Modifications', 
                        'Proteins': 'd_Proteins', 
                        'Gene Names': 'd_Gene Names',
                        'Protein Names': 'd_Protein Names', 
                        'Charge': 'd_Charge', 
                        'm/z': 'd_m/z',
                        'Mass': 'd_Mass', 
                        'Missed cleavages': 'd_missed_cleavages', 
                        'Length': 'd_length',
                        'Reverse': 'd_Reverse'}
    mg1.rename(columns=replacement_dict, inplace=True)
    grpdf = mg1.groupby('clusterID', as_index=False)[
        [modseq] + list(replacement_dict.values())].agg(
        lambda x: get_main_object(set(x)))
    grpdf.dropna(subset=[modseq], inplace=True)
    grpdf.drop(columns=[modseq], inplace=True)

    replacement_dict_reverse = { v : k for k, v in replacement_dict.items() }
    grpdf.rename(columns=replacement_dict_reverse,
                 inplace=True)
    grpdf[ident] = 't'
    mg1 = pd.merge(left=mg1, right=grpdf, on=['clusterID', ident], how='left')
    mg1.loc[
        mg1[ident] == 't', replacement_dict.values()] = mg1.loc[
        mg1[ident] == 't', replacement_dict.keys()].to_numpy()
    mg1.drop(columns=replacement_dict.keys(),
             inplace=True)
    mg1.rename(columns=replacement_dict_reverse, inplace=True)
    return mg1


def generate_summary_file(msmsscansdf, msmsdf, summarytxt, clusterfile):
    """
    Merges msmsscans.txt, msms.txt, and MaRaCluster clusters to generate summary file
    :param msmsscansdf:
    :param msmsdf:
    :param summarytxt:
    :param clusterfile:
    :return:
    """
    summary = merge_with_msmsscanstxt(clusterfile, msmsscansdf)
    summary = merge_with_summarytxt(summary, summarytxt)
    summary = merge_with_msmstxt(summary, msmsdf)
    return summary


def substitute_modifications(peps):
    """
    Substitutes p[STY] annotation with 1, 2, 3 to prevent miscounting of phosphoresidues
    :param peps:
    :return:
    """
    if type(peps) == str:
        rx = re.compile(r'pS')
        if re.search(rx, peps):
            peps = rx.sub('1', peps)
        rx = re.compile(r'pT')
        if re.search(rx, peps):
            peps = rx.sub('2', peps)
        rx = re.compile(r'pY')
        if re.search(rx, peps):
            peps = rx.sub('3', peps)
        return peps
    elif type(peps) == list:
        newlist = []
        for item in peps:
            item = substitute_modifications(item)
            newlist.append(item)
        return newlist
    elif peps != peps:
        return np.nan
    else:
        raise ValueError('String, string list, or np.NaN required for modified sequences')


def flag_ambiguous_clusters(sumdf, rawseq='Sequence', modseq='Modified sequence'):
    """
    adds MIC flags to all scans in MICs
    :param sumdf: Summary dataframe, merged from cleaned msms.txt and MaRaCluster clusters.tsv file
    :param rawseq: Column name of the raw sequence column
    :param modseq: Column name of the modified sequence column
    :return:
    """
    # Store number of unique modified sequences per cluster in group_dataframe
    group_dataframe = sumdf.groupby('clusterID')[modseq].nunique(dropna=True).reset_index(name='seqcount')

    # Filter for clusters with more than one modified sequence and flag with column mod_ambiguous
    group_dataframe = group_dataframe[group_dataframe['seqcount'] >= 2]
    group_dataframe.drop(columns=['seqcount'], axis=1, inplace=True)
    group_dataframe['mod_ambiguous'] = 1

    # merge with initial frame
    merged_dataframe = pd.merge(left=sumdf, right=group_dataframe, on='clusterID', how='left')

    # Store number of unique raw sequences per cluster in group_dataframe
    group_dataframe = sumdf.groupby('clusterID')[rawseq].nunique(dropna=True).reset_index(name='seqcount')

    # Filter for clusters with more than one raw sequence and flag with column raw_ambiguous
    group_dataframe = group_dataframe[group_dataframe['seqcount'] >= 2]
    group_dataframe.drop(columns=['seqcount'], axis=1, inplace=True)
    group_dataframe['raw_ambiguous'] = 1

    # merge with initial frame
    merged_dataframe = pd.merge(left=merged_dataframe, right=group_dataframe, on='clusterID', how='left')
    return merged_dataframe


def get_main_object(input_list):
    """
    Returns nan if no unique sequence is found in inputlist, or the sequence if its unique
    :param input_list: List of peptide sequences
    :return: Returns peptide sequence if it is the only one in the input list, otherwise returns np.NaN
    """
    x = remove_nan_values(input_list)
    if len(x) == 1:
        return x[0]
    else:
        return np.nan


def count_phos(modseq):
    if type(modseq) == str:
        return modseq.count('p')
    return np.nan


def count_clustering_parameters(summary, rawtrans=False):
    """
    Counts MICs, tIDs, dIDs, optionally IDs with lost phospho location, and clusters. Requires flagged MICs and tIDs.
    :param summary: summary_extended input dataframe
    :param rawtrans: Flag for added lost phospho localization counting
    :return: Dictionary of counted values
    """
    scans = len(summary)
    dids = len(summary[summary['identification'] == 'd'])
    tids = len(summary[summary['identification'] == 't'])
    lostphos = None
    if rawtrans:
        lostphos = len(
            summary[
                (summary['identification'] == 't') & (summary['Modified sequence'] != summary['Modified sequence'])])
    totclus = max(summary['clusterID'])
    mulclus = max(summary[summary['clusterID'].duplicated(keep=False)]['clusterID'])
    pho_mics = summary[summary['mod_ambiguous'] == 1]['clusterID'].nunique(dropna=True)
    raw_mics = summary[summary['raw_ambiguous'] == 1]['clusterID'].nunique(dropna=True)
    logger.info(f'Scans: {str(scans)}')
    logger.info(f'MaxQuant IDs: {str(dids)}')
    logger.info(f'SIMSI IDs: {str(tids)}')
    if rawtrans:
        logger.info(f'Isomeric SIMSI IDs: {str(lostphos)}')
    logger.info(f'All clusters: {str(totclus)}')
    logger.info(f'Clusters > 1: {str(mulclus)}')
    logger.info(f'PTM-isomeric clusters : {str(pho_mics)}')
    logger.info(f'Ambiguous clusters: {str(raw_mics)}')
    if rawtrans:
        return {'scans': scans, 'dids': dids, 'tids': tids, 'lostphos': lostphos, 'totclus': totclus,
                'mulclus': mulclus, 'pho_mics': pho_mics, 'raw_mics': raw_mics}
    else:
        return {'scans': scans, 'dids': dids, 'tids': tids, 'totclus': totclus, 'mulclus': mulclus,
                'pho_mics': pho_mics, 'raw_mics': raw_mics}


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


def assign_evidence_type(summary: pd.DataFrame, type_column_name: str = 'new_type'):   
    '''
    assign the updated Type column by checking if retention time of MS2 scan is within its evidence precursor rt window
    '''
    summary[type_column_name] = np.NaN
    msms_in_rt_window = ((summary['Retention time'] >= summary['Calibrated retention time start'] - summary['Retention time calibration']) & (
                          summary['Retention time'] <= summary['Calibrated retention time finish'] - summary['Retention time calibration']))
    summary.loc[msms_in_rt_window, type_column_name] = 'MULTI-MSMS'
    summary.loc[~msms_in_rt_window, type_column_name] = 'MSMS'
    
    # clear out the retention time columns for MS2 scans without precursor
    del_columns = ['Type', 'Calibrated retention time', 'Calibrated retention time start',
                   'Calibrated retention time finish']
    summary.loc[~msms_in_rt_window, del_columns] = np.NaN
    return summary


def assign_missing_precursors(summary: pd.DataFrame, allpeptides: pd.DataFrame):
    '''
    find precursors in allPeptides.txt for transferred MS2 scans in runs where the (peptide, charge) combination was previously not identified
    currently very slow (~500 MS2 scans / second). The apply function could be parallelized with the multiprocessing or Dask modules
    '''
    group_key = ['Raw file', 'Charge']
    allpeptides_grouped = { group : df_grouped for group, df_grouped in allpeptides.groupby(group_key) }
    
    missing_precursor = (summary['new_type'] == 'MSMS')
    summary.loc[missing_precursor, ['new_type', 'Intensity']] = summary.loc[missing_precursor].apply(
            lambda x : match_precursor_grouped(x, allpeptides_grouped, group_key), axis=1, result_type='expand'
        ).to_numpy() # need to convert to numpy array, see: https://stackoverflow.com/questions/69954697/
    return summary


def get_ppm_diff(mz1, mz2):
    return np.abs(mz1 - mz2)/mz1*1e6


def match_precursor_grouped(msms_scan: pd.Series, allpeptides_grouped: pd.core.groupby.DataFrameGroupBy, group_key: List[str]):
    if tuple(msms_scan[group_key]) not in allpeptides_grouped:
        return ['MSMS', np.NaN]
    allpeptides = allpeptides_grouped[tuple(msms_scan[group_key])]
    return match_precursor(msms_scan, allpeptides)


def match_precursor(msms_scan: pd.Series, allpeptides: pd.DataFrame, ppm_tol: float = 20.0):
    precursors = allpeptides[(get_ppm_diff(allpeptides['m/z'], msms_scan['m/z']) < ppm_tol) & \
            (allpeptides['Min scan number'] <= msms_scan['MS scan number']) & \
            (allpeptides['Max scan number'] >= msms_scan['MS scan number'])]
    if len(precursors.index) == 0:
        return ['MSMS', np.NaN]
    return ['MULTI-MSMS', precursors['Intensity'].values[0]]


def remove_duplicate_msms(summary: pd.DataFrame):
    '''
    remove duplicate MS2 scans due to (peptide, raw_file, charge) combinations matching multiple precursors
    prioritize the entry that matched a precursor (MULTI-MSMS) or has the lower retention time
    '''
    summary = summary.sort_values(by=['new_type', 'Calibrated retention time'], ascending=[False, True])
    summary = summary.drop_duplicates(subset=['summary_ID'], keep='first')
    return summary


def fill_missing_evidence_ids(summary):
    '''
    fill the evidence_ID column of MS2 scans without a precursor
    '''
    missing_evidence_id = summary['evidence_ID'].isna()
    
    start_id = int(summary['evidence_ID'].max() + 1)
    end_id = start_id + missing_evidence_id.sum()
    
    summary.loc[missing_evidence_id, 'evidence_ID'] = range(start_id, end_id)
    summary['evidence_ID'] = summary['evidence_ID'].astype(int)
    return summary


def assign_evidence_feature(summary: pd.DataFrame, evidence: pd.DataFrame, allpeptides: pd.DataFrame):
    # store number of rows of summary dataframe to check if we have the same number after merging
    summary_length_before_processing = len(summary.index)

    summary = merge_summary_with_evidence(summary, evidence)
    summary = assign_evidence_type(summary)
    summary = remove_duplicate_msms(summary)
    summary = assign_missing_precursors(summary, allpeptides)
    summary = fill_missing_evidence_ids(summary)
    
    if not len(summary) == summary_length_before_processing:
        raise ValueError(
            f'Number of summary entries changed during evidence feature assembly!'
            f'\n{summary_length_before_processing} before assembly,\n{len(summary)} after assembly.'
        )
    
    if not summary['summary_ID'].nunique() == len(summary):
        raise ValueError('Some summary_IDs were detected multiple times!')
    
    # TODO: check if all MULTI-MSMS entries are allocated to the correct PSMs
    
    summary = summary.sort_values(by=['Sequence', 'Modified sequence'])
    summary = summary.drop(columns=['Type'])
    return summary


def calculate_evidence_columns(summary, tmt):
    # replacing zeros with NaNs to count later
    reps = ['Reporter intensity 1', 'Reporter intensity 2', 'Reporter intensity 3', 'Reporter intensity 4',
            'Reporter intensity 5', 'Reporter intensity 6', 'Reporter intensity 7', 'Reporter intensity 8',
            'Reporter intensity 9', 'Reporter intensity 10', 'Reporter intensity corrected 1',
            'Reporter intensity corrected 2', 'Reporter intensity corrected 3', 'Reporter intensity corrected 4',
            'Reporter intensity corrected 5', 'Reporter intensity corrected 6', 'Reporter intensity corrected 7',
            'Reporter intensity corrected 8', 'Reporter intensity corrected 9', 'Reporter intensity corrected 10']
    if tmt == 11:
        reps.extend(['Reporter intensity 11', 'Reporter intensity corrected 11'])
    summary[reps] = summary[reps].replace({0: np.nan})
    summary = summary.sort_values(by=['Sequence', 'Modified sequence', 'Raw file', 'Charge']).reset_index(drop=True)

    summary_grouped = summary.groupby('evidence_ID')
    # Column generation
    # TODO: Sort for keeping best hit at top; best hit definition needed
    # TODO: Check every column; what is needed and is something missing?

    def csv_list(x):
        return ";".join(map(str, x))
    
    def csv_list_unique(x):
        return ";".join(map(str, list(dict.fromkeys(x))))

    evidence = summary_grouped.agg(
        **{
            'Sequence': pd.NamedAgg(column='Sequence', aggfunc='first'),
            'Length': pd.NamedAgg(column='Length', aggfunc='first'),
            'Modifications': pd.NamedAgg(column='Modifications', aggfunc='first'),
            'Modified sequence': pd.NamedAgg(column='Modified sequence', aggfunc='first'),
            'Missed cleavages': pd.NamedAgg(column='Missed cleavages', aggfunc='first'),
            'Proteins': pd.NamedAgg(column='Proteins', aggfunc=csv_list_unique), # from msms.txt
            'Leading proteins': pd.NamedAgg(column='Leading proteins', aggfunc=csv_list_unique), # from evidence.txt, NaN if scan not matched to precursor in evidence.txt
            'Gene Names': pd.NamedAgg(column='Gene Names', aggfunc='first'), # from msms.txt
            'Protein Names': pd.NamedAgg(column='Protein Names', aggfunc='first'),
            'Type': pd.NamedAgg(column='new_type', aggfunc='first'),
            'Raw file': pd.NamedAgg(column='Raw file', aggfunc='first'),
            'Fraction': pd.NamedAgg(column='Fraction', aggfunc='first'),
            'Experiment': pd.NamedAgg(column='Experiment', aggfunc='first'),
            'Charge': pd.NamedAgg(column='Charge', aggfunc='first'),
            'm/z': pd.NamedAgg(column='m/z', aggfunc='first'), # from msmsScans.txt, NaN if charge state in msmsScans is different than in the mzML
            'Mass': pd.NamedAgg(column='Mass', aggfunc='first'), # from msmsScans.txt, NaN if charge was unknown
            'Mass error [ppm]': pd.NamedAgg(column='Mass error [ppm]', aggfunc=min), # from msms.txt, NaN if Type=MSMS
            'Retention time': pd.NamedAgg(column='Retention time', aggfunc='first'), # from msmsScans.txt
            'PEP': pd.NamedAgg(column='PEP', aggfunc=min),
            'MS/MS count': pd.NamedAgg(column='Sequence', aggfunc='size'),
            'MS/MS all scan numbers': pd.NamedAgg(column='scanID', aggfunc=csv_list),
            'MS/MS scan number': pd.NamedAgg(column='scanID', aggfunc='first'),
            'Score': pd.NamedAgg(column='Score', aggfunc=max),
            'Delta score': pd.NamedAgg(column='Delta score', aggfunc=max),
            'Intensity': pd.NamedAgg(column='Intensity', aggfunc='sum'), # from evidence.txt, supplemented from allPeptides.txt     
            **{f'Reporter intensity corrected {i}': pd.NamedAgg(column=f'Reporter intensity corrected {i}', aggfunc='sum') for i in range(1, tmt+1)},
            **{f'Reporter intensity {i}': pd.NamedAgg(column=f'Reporter intensity {i}', aggfunc='sum') for i in range(1, tmt+1)},
            **{f'Reporter intensity count {i}': pd.NamedAgg(column=f'Reporter intensity {i}', aggfunc='count') for i in range(1, tmt+1)},
            'Reverse': pd.NamedAgg(column='Reverse', aggfunc='first'),
            'summary_ID': pd.NamedAgg(column='summary_ID', aggfunc=csv_list)
        })
    evidence['id'] = evidence.index # evidence_ID
    evidence = evidence.astype({'Length': 'int64', 'Missed cleavages': 'int64', 'Fraction': 'int64', 'Charge': 'int64', 'MS/MS scan number': 'int64'})
    evidence['Reverse'].fillna('', inplace=True)
    evidence = evidence.sort_values(by=['Sequence', 'Modified sequence', 'Raw file', 'Charge'])
    return evidence


def build_evidence(summary: pd.DataFrame, evidence: pd.DataFrame, allpeptides: pd.DataFrame, tmt: int):
    evidence = evidence[evidence['Type'] != 'MSMS']
    evidence = evidence.sort_values(by=['Sequence', 'Modified sequence', 'Raw file', 'Calibrated retention time start'])
    evidence.insert(len(evidence.columns), 'evidence_ID', range(len(evidence)))
    summary = assign_evidence_feature(summary, evidence, allpeptides)
    evidence = calculate_evidence_columns(summary, tmt)
    return evidence


def remove_unidentified_scans(summary):
    summary = summary.loc[~summary['identification'].isna()]
    summary.insert(0, 'summary_ID', range(len(summary)))
    return summary


if __name__ == '__main__':
    raise NotImplementedError('Do not run this script.')
