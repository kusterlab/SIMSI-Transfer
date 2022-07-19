import io

import pytest
import pandas as pd
import numpy as np

import simsi_transfer.transfer as transfer


pd.set_option('display.max_columns', None)

class TestTransfer:
    def test_transfer(self, summary_df):
        transferred_df = transfer.transfer(summary_df, ambiguity_decision=False)
        assert transferred_df.iloc[3]['Sequence'] == 'DSDSWDADAFSVEDPVRK'
        assert transferred_df.iloc[3]['Modified sequence'] == 'DSDSWDADAFSVEDPVRK.2.p2/p4/p11'

    def test_transfer_ambiguity_decision(self, summary_df):
        transferred_df = transfer.transfer(summary_df, ambiguity_decision=True)
        assert transferred_df.iloc[3]['Sequence'] == 'DSDSWDADAFSVEDPVRK'
        assert transferred_df.iloc[3]['Modified sequence'] == 'DS(Phospho (STY))DS(Phospho (STY))WDADAFSVEDPVRK'


def test_get_modified_sequence_annotation():
    test1 = ['SSS(Phospho (STY))PPPRK', 'SS(Phospho (STY))SPPPRK', 'SSS(Phospho (STY))PPPRK', 'SS(Phospho (STY))SPPPRK']
    test1_psp_format = ['SSsPPPRK', 'SsSPPPRK', 'SSsPPPRK', 'SsSPPPRK']
    assert transfer.get_modified_sequence_annotation(test1, test1_psp_format) == 'SSSPPPRK.1.p2/p3'


def TestCalculateAverageProbabilities():
    def test_calculate_average_probabilities(self):
        test1 = ['S(0.326)S(0.185)S(0.489)PPPRK', 'SS(1)SPPPRK', 'S(0.01)S(0.19)S(0.8)PPPRK', 'S(0.4)S(0.6)SPPPRK']
        assert transfer.calculate_average_probabilities(test1) == 'S(0.184)S(0.494)S(0.322)PPPRK'
    
    def test_calculate_average_probabilities_transfers(self):
        test2 = [np.nan, np.nan, np.nan, 'S(0.4)S(0.6)SPPPRK']
        assert transfer.calculate_average_probabilities(test2) == 'S(0.4)S(0.6)SPPPRK'
    
    def test_calculate_average_probabilities_nan(self):
        test3 = [np.nan]
        assert np.isnan(transfer.calculate_average_probabilities(test3))


def test_remove_probabilities_from_sequence():
    sequence_with_probabilities = 'AAAAAAAATMALAAPS(0.059)S(0.899)PT(0.176)PES(0.701)PT(0.163)MLT(0.001)K'
    assert transfer.remove_probabilities_from_sequence(sequence_with_probabilities) == 'AAAAAAAATMALAAPSSPTPESPTMLTK'


def test_remove_probabilities_from_sequence_integer():
    sequence_with_probabilities = 'AAAAAAAGDS(1)DS(1)WDADAFSVEDPVRK'
    assert transfer.remove_probabilities_from_sequence(sequence_with_probabilities) == 'AAAAAAAGDSDSWDADAFSVEDPVRK'


def test_remove_probabilities_from_sequence_no_mods():
    sequence_with_probabilities = 'AAAAAAAGDSDSWDADAFSVEDPVRK'
    assert transfer.remove_probabilities_from_sequence(sequence_with_probabilities) == 'AAAAAAAGDSDSWDADAFSVEDPVRK'


def test_get_mod_probabilities_dict():
    sequence_with_probabilities = 'AAAAAAAATMALAAPS(0.059)S(0.899)PT(0.176)PES(0.701)PT(0.163)MLT(0.001)K'
    assert transfer.get_mod_probabilities_dict(sequence_with_probabilities) == {16: 0.059, 17: 0.899, 19: 0.176, 22: 0.701, 24: 0.163, 27: 0.001}


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
    df['Modifications'] = ""
    df['Proteins'] = ""
    df['Gene Names'] = ""
    df['Protein Names'] = ""
    df['Charge'] = 2
    df['m/z'] = 500.0
    df['Mass'] = 1000.0
    df['Missed cleavages'] = 0
    df['Length'] = 1
    df['Reverse'] = ""
    return df

