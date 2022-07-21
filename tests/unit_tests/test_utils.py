import pandas as pd
import numpy as np

import simsi_transfer.utils.utils as utils


def test_apply_and_flatten():
    raw_folders = ['folder_1', 'folder_2', 'folder_3']
    def list_files(folder):
        return [f'{folder}/file_1.txt', f'{folder}/file_1.txt', f'{folder}/file_2.txt']
    
    raw_files = utils.apply_and_flatten(raw_folders, list_files)
    assert raw_files == ['folder_1/file_1.txt', 'folder_1/file_2.txt', 'folder_2/file_1.txt', 'folder_2/file_2.txt', 'folder_3/file_1.txt', 'folder_3/file_2.txt']


def test_csv_list_unique():
    x = pd.Series(['ABC','BCD','BCD', np.nan], name='Proteins')
    assert utils.csv_list_unique(x) == 'ABC;BCD'


def test_csv_list_unique_multiple_isoforms():
    x = pd.Series(['ABC;BCD','BCD;ABC','BCD'], name='Proteins')
    assert utils.csv_list_unique(x) == 'ABC;BCD'
