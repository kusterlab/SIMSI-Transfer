import re
import logging

import pandas as pd
import numpy as np

from .merging_functions import merge_with_msmsscanstxt, merge_with_msmstxt, merge_with_summarytxt


logger = logging.getLogger(__name__)


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
        name='cluster_modified_sequences')

    grpdf[['representative_modified_sequence', 'representative_raw_sequence']] = grpdf[
        'cluster_modified_sequences'].apply(get_modified_and_raw_sequence)

    grpdf.loc[grpdf['representative_modified_sequence'].notna(), ident] = 't'
    mg1 = pd.merge(left=sumdf, right=grpdf, on=['clusterID'], how='left')
    mg1.loc[mg1[modseq] == mg1[modseq], ident] = 'd'
    mg1.loc[mg1[modseq] != mg1[modseq], modseq] = mg1['representative_modified_sequence']
    mg1.loc[mg1[modseq] != mg1[modseq], rawseq] = mg1['representative_raw_sequence']
    mg1.drop(columns=['representative_modified_sequence', 'representative_raw_sequence'], axis=1, inplace=True)

    replacement_dict = {'Modifications': 'd_Modifications',
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
    grpdf = mg1.groupby('clusterID', as_index=False)[[rawseq] + list(replacement_dict.values())].agg(
        lambda x: get_main_object(set(x)))
    grpdf.dropna(subset=[rawseq], inplace=True)
    grpdf.drop(columns=[rawseq], inplace=True)

    replacement_dict_reverse = {v: k for k, v in replacement_dict.items()}
    grpdf.rename(columns=replacement_dict_reverse, inplace=True)
    grpdf[ident] = 't'
    mg1 = pd.merge(left=mg1, right=grpdf, on=['clusterID', ident], how='left')
    mg1.loc[
        mg1[ident] == 't', replacement_dict.values()] = mg1.loc[
        mg1[ident] == 't', replacement_dict.keys()].to_numpy()
    mg1.drop(columns=replacement_dict.keys(), inplace=True)
    mg1.rename(columns=replacement_dict_reverse, inplace=True)
    return mg1


def get_modified_and_raw_sequence(sequence_array):
    if len(sequence_array) == 0:
        return pd.Series([np.NaN, np.NaN])

    raw_sequence = remove_modifications(sequence_array[0])
    if len(sequence_array) == 1:
        return pd.Series([sequence_array[0], raw_sequence])

    sequence_array_lower = [re.sub(re.compile(r'([STY])\(Phospho \(STY\)\)'), lambda pat: pat.group(1).lower(), x) for x in sequence_array]
    sequence_set_upper = {x.upper() for x in sequence_array_lower}
    if len(sequence_set_upper) != 1:
        return pd.Series([np.NaN, np.NaN])

    modified_sequence = generate_modified_sequence_annotation(sequence_array, sequence_array_lower)
    return pd.Series([modified_sequence, raw_sequence])


def generate_modified_sequence_annotation(sequence_array, sequence_array_lower):
    num_mods = len(re.findall(r'[sty]', sequence_array_lower[0]))
    phospho_positions = {m.start() + 1 for x in sequence_array_lower for m in re.finditer(r'[sty]', x)}
    phospho_positions = sorted(list(phospho_positions))
    phopsho_positions_string = "/".join(map(lambda x: f'p{x}', phospho_positions))
    mod_sequence_without_phospho = remove_modifications(sequence_array[0], remove_phospho_only=True)
    return f"{mod_sequence_without_phospho}.{num_mods}.{phopsho_positions_string}"


def remove_modifications(mod_sequence, remove_phospho_only=False):
    raw_sequence = re.sub(re.compile(r'([STY])\(Phospho \(STY\)\)'), '', mod_sequence)
    if remove_phospho_only:
        return raw_sequence
    raw_sequence = raw_sequence.replace('(Acetyl (Protein N-term))', '')
    raw_sequence = raw_sequence.replace('(Oxidation (M))', '')
    raw_sequence = raw_sequence.replace('_', '')
    return raw_sequence


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


def remove_nan_values(input_list):
    """
    Eliminates NaNs from sets or lists, returns list of all other values
    :param input_list: List (or set) of elements
    :return: List with removed np.NaN values
    """
    return [v for v in input_list if pd.notnull(v)]


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