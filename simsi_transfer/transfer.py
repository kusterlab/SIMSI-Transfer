import collections
import functools
import operator
import re
import logging
from typing import List, Callable, Dict

import pandas as pd
import numpy as np

from .utils import utils

logger = logging.getLogger(__name__)

PHOSPHO_REGEX = re.compile(r'([STY])\(Phospho \(STY\)\)')
PROBABILITY_REGEX = re.compile(r'\((\d(?:\.?\d+)?)\)')


def transfer(summary_df, max_pep=False, mask=False, ambiguity_decision='majority'):
    """
    Main function for transfers by clustering. Transfers identifications for merged dataframe and adds a column for
    identification type. Transferred columns are Sequence, Modified sequence, Proteins, Gene names, Protein Names,
    Charge, m/z, and Mass.
    :param summary_df: Summary dataframe, merged from cleaned msms.txt and MaRaCluster clusters.tsv file
    :param max_pep: Maximum PEP of a PSM to be considered for transfers
    :param mask: if false, uses "identification" column from summary frame. If set to a value, transfer uses the
    'identification_{mask}' column, needed for the masking analysis
    :param ambiguity_decision: if 'all', returns raw sequences with potential phospho positions when encountering isomer
    clusters; if 'majority', decides for one sequence by majority vote
    :return: DataFrame with transferred identifications resembling MaxQuant msmsScans.txt
    """
    if mask:
        identification_column = f'identification_{mask}'
    else:
        identification_column = 'identification'

    if ambiguity_decision == 'majority':
        # TODO: Generate modified sequence from probability string rather than taking it from the cluster
        mod_seq_func = lambda s: get_consensus_modified_sequence(s, get_most_common_sequence)
    elif ambiguity_decision == 'all':
        mod_seq_func = get_consensus_modified_sequence
    else:
        raise ValueError("The parameter 'ambiguity_decision' has to be set on 'all' or 'majority'!")

    agg_funcs = {'Sequence': utils.get_unique_else_nan,
                 'Modifications': utils.get_unique_else_nan,
                 'Modified sequence': mod_seq_func,
                 'Phospho (STY) Probabilities': calculate_average_probabilities,
                 'Proteins': utils.csv_list_unique,
                 'Gene Names': utils.csv_list_unique,
                 'Protein Names': utils.csv_list_unique,
                 'Charge': utils.get_unique_else_nan,
                 'm/z': 'mean',
                 'Mass': 'mean',
                 'Missed cleavages': utils.get_unique_else_nan,
                 'Length': utils.get_unique_else_nan,
                 'PEP': 'max',
                 'Reverse': utils.get_unique_else_nan}

    identified_scans = summary_df['Modified sequence'].notna()
    pep_filtered = pd.Series(np.ones_like(summary_df.index), dtype='bool')
    if max_pep:
        pep_filtered = summary_df['PEP'].astype(float) <= max_pep / 100
    cluster_info_df = summary_df[identified_scans & pep_filtered].groupby('clusterID', as_index=False).agg(agg_funcs)
    # Mark all clusters with a unique identification as transferred ('t').
    # Identifications by MQ will overwrite this column as direct identification ('d') a few lines below.
    cluster_info_df[identification_column] = np.where(cluster_info_df['Modified sequence'].notna(), 't', None)

    replacement_dict = {k: "grouped_" + k for k in agg_funcs.keys()}
    cluster_info_df.rename(columns=replacement_dict, inplace=True)

    summary_df = pd.merge(left=summary_df, right=cluster_info_df, on=['clusterID'], how='left')
    summary_df.loc[identified_scans, identification_column] = 'd'
    # Now we have 'd' in every ID that came from MaxQuant and 't' in every ID that came from the clustering
    # And now we copy the contents of the 'grouped_...' columns into the original columns for every row where we have a 't'
    summary_df.loc[
        # Aren't all of these cells by definition NaNs? Can we make use of that?
        # Like 'Merge, and if you find two columns with the same name always use the value that is not NaN!'
        # Then we could do this directly in the merge and wouldn't have to go over the df again here
        # Try pd.DataFrame.combine(). For this the dfs would need to be in the same shape though.
        summary_df[identification_column] == 't', replacement_dict.keys()] = summary_df.loc[
        summary_df[identification_column] == 't', replacement_dict.values()].to_numpy()
    # The 'grouped_...' columns have made themselves redundant
    summary_df.drop(columns=replacement_dict.values(), inplace=True)
    return summary_df


def transform_phospho_psp_format(sequences: List[str]) -> List[str]:
    """Replaces phosphorylation modification by lower case letter (PhosphositePlus convention), 
    e.g. 'AAAAAAAGDS(Phospho (STY))DS(Phospho (STY))WDADAFSVEDPVRK' => 'AAAAAAAGDsDsWDADAFSVEDPVRK'
    """
    sequences = [re.sub('_', '', x) for x in sequences]
    return [re.sub(PHOSPHO_REGEX, lambda pat: pat.group(1).lower(), x) for x in sequences]


def remove_probabilities(sequences_with_probabilities: List[str]) -> List[str]:
    return [remove_probabilities_from_sequence(s) for s in sequences_with_probabilities]


def remove_probabilities_from_sequence(sequence: str) -> str:
    """remove probability strings noted in parentheses, 
    e.g. 'AAAAAAAGDS(0.988)DS(0.012)WDADAFSVEDPVRK' => 'AAAAAAAGDSDSWDADAFSVEDPVRK'
    """
    return re.sub(PROBABILITY_REGEX, '', sequence)


def check_ambiguity(sequences: List[str],
                    transform_sequence: Callable[[List[str]], List[str]] = transform_phospho_psp_format):
    sequences = utils.remove_nan_values(set(sequences))
    if len(sequences) == 0:
        return np.nan, np.nan
    if len(sequences) == 1:
        return sequences[0], np.nan

    transformed_sequences = transform_sequence(sequences)
    unique_sequences = set(map(str.upper, transformed_sequences))
    if len(unique_sequences) != 1:
        return np.NaN, np.nan
    return sequences, transformed_sequences


def is_unambiguous(sequences):
    return (type(sequences) == float and np.isnan(sequences)) or type(sequences) == str


def get_mod_probabilities_dict(probstring: str) -> Dict[int, float]:
    probsplit = re.split(PROBABILITY_REGEX, probstring)
    start = 0
    probdict = dict()
    for amino_acids, prob in zip(probsplit[0::2], probsplit[1::2]):
        start += len(amino_acids)
        probdict[start] = float(prob)
    return probdict


def sum_dictionaries(dictionary):
    return dict(functools.reduce(operator.add, map(collections.Counter, dictionary)))


def average_dictionaries(dictionary):
    summed_dictionary = sum_dictionaries(dictionary)

    num_dictionaries = len(dictionary)

    def get_average_and_round(x):
        return round(x / num_dictionaries, 3)

    return {pos: get_average_and_round(summed_prob) for pos, summed_prob in summed_dictionary.items()}


def add_probabilities_to_sequence(sequence, probability_dict):
    for position, probability in reversed(sorted(list(probability_dict.items()))):
        sequence = sequence[:position] + f'({probability})' + sequence[position:]
    return sequence


def calculate_average_probabilities(mod_probability_sequences):
    mod_probability_sequences, _ = check_ambiguity(mod_probability_sequences, remove_probabilities)
    if is_unambiguous(mod_probability_sequences):
        return mod_probability_sequences

    mod_probabilities = [get_mod_probabilities_dict(p) for p in mod_probability_sequences]
    average_mod_probabilities = average_dictionaries(mod_probabilities)

    sequence = remove_probabilities_from_sequence(mod_probability_sequences[0])
    sequence_with_probabilities = add_probabilities_to_sequence(sequence, average_mod_probabilities)
    # # this does not work for multiple phosphorylations in one sequence
    # if not 0.95 < sum(resdict.values()) < 1.05:
    #     raise ValueError(f'Probability sum deviates from expected value: {resdict}')
    return sequence_with_probabilities


def get_modified_sequence_annotation(sequences, sequences_psp_format):
    # TODO: Make ambiguous localization handling work for acetylation (and p-ac-combinations) as well
    # TODO: If underscores are handed over in the modified seq, the modification positions are incorrect in the end!
    num_mods = len(re.findall(r'[sty]', sequences_psp_format[0]))
    phospho_positions = {m.start() + 1 for x in sequences_psp_format for m in re.finditer(r'[sty]', x)}
    phospho_positions = sorted(list(phospho_positions))
    phopsho_positions_string = "/".join(map(lambda x: f'p{x}', phospho_positions))
    mod_sequence_without_phospho = utils.remove_modifications(sequences[0], remove_phospho_only=True)
    return f"{mod_sequence_without_phospho}.{num_mods}.{phopsho_positions_string}"


def get_most_common_sequence(sequences, sequences_psp_format):
    sequences = utils.remove_nan_values(sequences)
    return collections.Counter(sequences).most_common(1)[0][0]


def get_consensus_modified_sequence(sequences, consensus_function=get_modified_sequence_annotation):
    unique_sequences, sequences_psp_format = check_ambiguity(sequences)
    if is_unambiguous(unique_sequences):
        return unique_sequences
    return consensus_function(sequences, sequences_psp_format)


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
