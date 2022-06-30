import collections
import functools
import operator
import re
import logging

import pandas as pd
import numpy as np

from .merging_functions import merge_with_msmsscanstxt, merge_with_msmstxt, merge_with_summarytxt


logger = logging.getLogger(__name__)

PHOSPHO_REGEX = re.compile(r'([STY])\(Phospho \(STY\)\)')


def transfer(summary_df, mask=False, ambiguity_decision=False):
    """
    Main function for transfers by clustering. Transfers identifications for merged dataframe and adds a column for
    identification type. Transferred columns are Sequence, Modified sequence, Proteins, Gene names, Protein Names,
    Charge, m/z, and Mass.
    :param summary_df: Summary dataframe, merged from cleaned msms.txt and MaRaCluster clusters.tsv file
    :param mask: if false, uses "identification" column from summary frame. If set to a value, transfer uses the
    'identification_{mask}' column, needed for the masking analysis
    :param ambiguity_decision: if False, returns raw sequences with potential phospho positions when encountering isomer
    clusters; if True, decides for one sequence (majority vote)
    :return: DataFrame with transferred identifications resembling MaxQuant msmsScans.txt
    """
    if mask:
        identification_column = f'identification_{mask}'
    else:
        identification_column = 'identification'

    identified_scans = summary_df['Modified sequence'].notna()

    # TODO: Generate modified sequence from probability string rather than taking it from the cluster
    agg_funcs = {'Sequence': get_unique_else_nan,
                 'Modifications': get_unique_else_nan,
                 'Phospho (STY) Probabilities': calculate_average_probabilities,
                 'Proteins': get_unique_else_nan,
                 'Gene Names': get_unique_else_nan,
                 'Protein Names': get_unique_else_nan,
                 'Charge': get_unique_else_nan,
                 'm/z': 'mean',
                 'Mass': 'mean',
                 'Missed cleavages': get_unique_else_nan,
                 'Length': get_unique_else_nan,
                 'Reverse': get_unique_else_nan}
    if ambiguity_decision:
        agg_funcs['Modified sequence'] = get_modified_sequence_decision
    else:
        agg_funcs['Modified sequence'] = get_consensus_modified_sequence

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



def check_ambiguity(sequences, checktype='sequences'):
    sequences = remove_nan_values(set(sequences))
    if len(sequences) == 0:
        return np.nan, np.nan
    if len(sequences) == 1:
        return sequences[0], np.nan

    sequences_lower = [re.sub(PHOSPHO_REGEX, lambda pat: pat.group(1).lower(), x) for x in sequences]
    if checktype == 'sequences':
        sequence_set_upper = {x.upper() for x in sequences_lower}
        if len(sequence_set_upper) != 1:
            return np.NaN, np.nan
    elif checktype == 'probabilities':
        pattern = re.compile(r'\((\d(?:\.\d+)?)\)')
        pureseq = {re.sub(pattern, '', sequence) for sequence in sequences}
        if len(pureseq) != 1:
            return np.NaN, np.nan
    return sequences, sequences_lower


def calculate_average_probabilities(initial_phospho_probabilities):
    phospho_probabilities, sequences_lower = check_ambiguity(initial_phospho_probabilities, 'probabilities')
    # logger.info(phospho_probabilities)
    if type(phospho_probabilities) == float:
        if np.isnan(phospho_probabilities):
            return np.nan
    elif type(phospho_probabilities) == str:
        return phospho_probabilities

    phospho_probabilities = remove_nan_values(initial_phospho_probabilities)

    pattern = re.compile(r'\((\d(?:\.?\d+)?)\)')
    pureseq = re.sub(pattern, '', phospho_probabilities[0])
    tots = len(phospho_probabilities)
    dictlist = []
    for probstring in phospho_probabilities:
        probsplit = re.split(pattern, probstring)
        start = 0
        probdict = dict()
        for listpos in range(len(probsplit)):
            if (listpos % 2 == 0) & (listpos + 1 != len(probsplit)):
                start = start + len(probsplit[listpos])
                probdict[start] = float(probsplit[listpos + 1])
        dictlist.append(probdict)
    resdict = dict(functools.reduce(operator.add, map(collections.Counter, dictlist)))
    resdict = dict(sorted({k: round(v / tots, 3) for k, v in resdict.items()}.items()))
    for key in reversed(list(resdict.keys())):
        pureseq = pureseq[:key] + f'({resdict[key]})' + pureseq[key:]
    # # this does not work for multiple phosphorylations in one sequence
    # if not 0.95 < sum(resdict.values()) < 1.05:
    #     raise ValueError(f'Probability sum deviates from expected value: {resdict}')
    return pureseq


def get_consensus_modified_sequence(sequences):
    sequences, sequences_lower = check_ambiguity(sequences)
    if type(sequences) == float:
        if np.isnan(sequences):
            return np.nan
    elif type(sequences) == str:
        return sequences

    modified_sequence = generate_modified_sequence_annotation(sequences, sequences_lower)
    return modified_sequence


def get_modified_sequence_decision(initial_sequences):
    sequences, sequences_lower = check_ambiguity(initial_sequences)
    if type(sequences) == float:
        if np.isnan(sequences):
            return np.nan
    elif type(sequences) == str:
        return sequences

    sequences = remove_nan_values(initial_sequences)
    return collections.Counter(sequences).most_common(1)[0][0]


def generate_modified_sequence_annotation(sequences, sequences_lower):
    # TODO: Make ambiguous localization handling work for acetylation (and p-ac-combinations) as well
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

    # remove mod_ambiguous flags if cluster is raw_ambiguous
    merged_dataframe.loc[merged_dataframe['raw_ambiguous'] == 1, 'mod_ambiguous'] = np.nan
    return merged_dataframe
