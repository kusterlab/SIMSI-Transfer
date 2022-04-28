import re
import logging

import pandas as pd
import numpy as np

from .merging_functions import merge_with_msmsscanstxt, merge_with_msmstxt, merge_with_summarytxt


logger = logging.getLogger(__name__)

PHOSPHO_REGEX = re.compile(r'([STY])\(Phospho \(STY\)\)')


def transfer(summary_df, mask=False):
    """
    Main function for transfers by clustering. Transfers identifications for merged dataframe and adds a column for
    identification type. Transferred columns are Sequence, Modified sequence, Proteins, Gene names, Protein Names,
    Charge, m/z, and Mass.
    :param summary_df: Summary dataframe, merged from cleaned msms.txt and MaRaCluster clusters.tsv file
    :param mask: if false, uses "identification" column from summary frame. If set to a value, transfer uses the
    'identification_{mask}' column, needed for the masking analysis
    :return: DataFrame with transferred identifications resembling MaxQuant msmsScans.txt
    """
    if mask:
        identification_column = f'identification_{mask}'
    else:
        identification_column = 'identification'

    identified_scans = summary_df['Modified sequence'].notna()

    agg_funcs = {'Sequence': get_unique_else_nan,
                 'Modified sequence': get_consensus_modified_sequence,
                 'Modifications': get_unique_else_nan,
                 'Proteins': get_unique_else_nan,
                 'Gene Names': get_unique_else_nan,
                 'Protein Names': get_unique_else_nan,
                 'Charge': get_unique_else_nan,
                 'm/z': 'mean',
                 'Mass': 'mean',
                 'Missed cleavages': get_unique_else_nan,
                 'Length': get_unique_else_nan,
                 'Reverse': get_unique_else_nan}

    cluster_info_df = summary_df[identified_scans].groupby('clusterID', as_index=False).agg(agg_funcs)
    
    # Mark all clusters with a unique identification as transferred ('t').
    # Identifications by MQ will overwrite this column as direct identification ('d') a few lines below.
    cluster_info_df[identification_column] = np.where(cluster_info_df['Modified sequence'].notna(), 't', None)
    
    replacement_dict = {k: "d_" + k for k in agg_funcs.keys()}
    cluster_info_df.rename(columns=replacement_dict, inplace=True)
    
    summary_df = pd.merge(left=summary_df, right=cluster_info_df, on=['clusterID'], how='left')
    summary_df.loc[identified_scans, identification_column] = 'd'
    summary_df.loc[
        summary_df[identification_column] == 't', replacement_dict.keys()] = summary_df.loc[
        summary_df[identification_column] == 't', replacement_dict.values()].to_numpy()
    summary_df.drop(columns=replacement_dict.values(), inplace=True)
    return summary_df


def get_consensus_modified_sequence(sequences):
    sequences = remove_nan_values(set(sequences))
    if len(sequences) == 0:
        return np.NaN
    
    if len(sequences) == 1:
        return sequences[0]

    sequences_lower = [re.sub(PHOSPHO_REGEX, lambda pat: pat.group(1).lower(), x) for x in sequences]
    sequence_set_upper = {x.upper() for x in sequences_lower}
    if len(sequence_set_upper) != 1:
        return np.NaN

    modified_sequence = generate_modified_sequence_annotation(sequences, sequences_lower)
    return modified_sequence


def generate_modified_sequence_annotation(sequences, sequences_lower):
    num_mods = len(re.findall(r'[sty]', sequences_lower[0]))
    phospho_positions = {m.start() + 1 for x in sequences_lower for m in re.finditer(r'[sty]', x)}
    phospho_positions = sorted(list(phospho_positions))
    phopsho_positions_string = "/".join(map(lambda x: f'p{x}', phospho_positions))
    mod_sequence_without_phospho = remove_modifications(sequences[0], remove_phospho_only=True)
    return f"{mod_sequence_without_phospho}.{num_mods}.{phopsho_positions_string}"


def remove_modifications(mod_sequence, remove_phospho_only=False):
    raw_sequence = re.sub(PHOSPHO_REGEX, '', mod_sequence)
    if remove_phospho_only:
        return raw_sequence
    raw_sequence = raw_sequence.replace('(Acetyl (Protein N-term))', '')
    raw_sequence = raw_sequence.replace('(Oxidation (M))', '')
    raw_sequence = raw_sequence.replace('_', '')
    return raw_sequence


def get_unique_else_nan(input_list):
    """
    Returns nan if no unique sequence is found in inputlist, or the sequence if its unique
    :param input_list: List of peptide sequences
    :return: Returns peptide sequence if it is the only one in the input list, otherwise returns np.NaN
    """
    x = remove_nan_values(set(input_list))
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
