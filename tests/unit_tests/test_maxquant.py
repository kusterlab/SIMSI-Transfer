import io

import pytest
import pandas as pd

import simsi_transfer.utils.utils as mq

pd.set_option('display.max_columns', None)


def test_process_and_concat(summary_df):
    mq_txt_folders = ['folder_1', 'folder_2', 'folder_3']
    def read_txt(folder):
        df = summary_df.copy()
        df['Batch'] = folder
        return df

    merged_df = mq.process_and_concat(mq_txt_folders, read_txt)
    assert len(merged_df.index) == 7*3


# Creating dataframes from strings: https://towardsdatascience.com/67b0c2b71e6a
@pytest.fixture
def summary_df():
    df_string = """clusterID; Sequence;               Modified sequence;                             Phospho (STY) Probabilities;
1; DSDSWDADAFSVEDPVRK; DS(Phospho (STY))DS(Phospho (STY))WDADAFSVEDPVRK;                                DS(1)DS(1)WDADAFSVEDPVRK;
1; DSDSWDADAFSVEDPVRK; DS(Phospho (STY))DS(Phospho (STY))WDADAFSVEDPVRK;                        DS(0.955)DS(0.045)WDADAFSVEDPVRK;
1; DSDSWDADAFSVEDPVRK; DS(Phospho (STY))DSWDADAFS(Phospho (STY))VEDPVRK;                        DS(0.765)DSWDADAFS(0.235)VEDPVRK;
1;                   ;                                                 ;                                                        ;
2;      SSPTPESPTMLTK;      SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK; S(0.059)S(0.899)PT(0.176)PES(0.701)PT(0.163)MLT(0.001)K;
2;                   ;                                                 ;                                                        ;
2;      SSPTPESPTMLTK;      SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK; S(0.005)S(0.205)PT(0.167)PES(0.526)PT(0.092)MLT(0.004)K;"""
    df = pd.read_csv(io.StringIO(df_string), delimiter=';', skipinitialspace=True)
    return df
