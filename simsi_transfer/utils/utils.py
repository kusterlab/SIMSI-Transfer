from typing import List, Callable, Any, Union, Set
from pathlib import Path
import logging

import pandas as pd
import numpy as np
from job_pool.tqdm_logger import TqdmToLogger
from tqdm import tqdm


logger = logging.getLogger(__name__)


def csv_list_unique(x: pd.Series) -> Union[str, float]:
    """Returns a unique semicolon-separated list from merging multiple semicolon-separated lists."""
    x = remove_nan_values(set(x))
    if len(x) == 1:
        return x[0]
    elif len(x) == 0:
        return np.nan

    def semicolon_split(x):
        return x.split(";")
    return ";".join(apply_and_flatten(x, semicolon_split))


def csv_unique(s: str) -> str:
    return ";".join(sorted(set([x for x in s.split(";") if len(x) > 0])))


def convert_to_path_list(s):
    return list(map(Path, s.tolist()))


def get_raw_files_and_correction_factor_paths(meta_input_df: pd.DataFrame):
    meta_input_exploded_df = meta_input_df[['raw_files', 'tmt_correction_file']].explode('raw_files')
    raw_files, correction_factor_paths = meta_input_exploded_df.transpose().values.tolist()
    
    def convert_to_paths(l: List[str]):
        return list(map(Path, l))

    raw_files = convert_to_paths(raw_files)
    correction_factor_paths = convert_to_paths(correction_factor_paths)
    
    return raw_files, correction_factor_paths


def apply_and_flatten(lst: List[Any], f: Callable[[Any],List[Any]]) -> List[Any]:
    """Applies a function that returns a list on each item of the list and flattens the result into a single list with unique items"""
    return sorted(list({f_x for x in lst for f_x in f(x)}))


def get_unique_else_nan(input_list: List[Any]) -> Any:
    """
    Returns nan if no unique sequence is found in inputlist, or the sequence if its unique
    :param input_list: List of peptide sequences
    :return: Returns peptide sequence if it is the only one in the input list, otherwise returns np.NaN
    """
    x = remove_nan_values(set(input_list))
    if len(x) == 1:
        return x[0]
    return np.nan


def remove_nan_values(input_list: Union[List[Any], Set[Any]]) -> List[Any]:
    """
    Eliminates NaNs from sets or lists, returns list of all other values
    :param input_list: List (or set) of elements
    :return: List with removed np.NaN values
    """
    return [v for v in input_list if pd.notnull(v)]


def remove_modifications(mod_sequence: str, remove_phospho_only: bool = False) -> str:
    raw_sequence = mod_sequence.replace('(Phospho (STY))', '')
    if remove_phospho_only:
        return raw_sequence
    raw_sequence = raw_sequence.replace('(Acetyl (Protein N-term))', '')
    raw_sequence = raw_sequence.replace('(Oxidation (M))', '')
    raw_sequence = raw_sequence.replace('_', '')
    return raw_sequence


def human_readable_size(size, decimal_places=2):
    for unit in ['B', 'KiB', 'MiB', 'GiB', 'TiB', 'PiB']:
        if size < 1024.0 or unit == 'PiB':
            break
        size /= 1024.0
    return f"{size:.{decimal_places}f} {unit}"


def get_dataframe_size(df: pd.DataFrame):
    return human_readable_size(df.memory_usage(deep=True).sum())


def process_and_concat(input_folders: List[Any], reading_function: Callable, **kwargs) -> pd.DataFrame:
    tqdm_out = TqdmToLogger(logger, level=logging.INFO)
    return pd.concat(
        [reading_function(f, **kwargs) for f in tqdm(input_folders, file=tqdm_out, mininterval=10)],
        axis=0,
    )