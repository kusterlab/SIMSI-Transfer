import io

import pytest
import pandas as pd
import pandas.util.testing as tm
import numpy as np

import simsi_transfer.evidence as evidence
import simsi_transfer.transfer as transfer


def test_assign_missing_precursors(summary_missing_precursors, allpeptides_missing_precursors):
    summary = evidence.assign_missing_precursors(summary_missing_precursors, allpeptides_missing_precursors)
    
    def check_cols(idx, vals):
        tm.assert_series_equal(summary.loc[idx, ['new_type', 'Intensity']].reset_index(drop=True), pd.Series(vals, name=idx))  
  
    # precursor found previously
    check_cols(0, ['MULTI-MSMS', 1238800.0])
    # precursor found
    check_cols(1, ['MULTI-MSMS', 2348800.0])
    # outside of scan range
    check_cols(2, ['MSMS', np.nan])
    # wrong file
    check_cols(3, ['MSMS', np.nan])
    # wrong charge
    check_cols(4, ['MSMS', np.nan])
    # missing file
    check_cols(5, ['MSMS', np.nan])

def test_generate_modified_sequence_annotation_unique_raw_sequence():
    assert (transfer.get_modified_and_raw_sequence(['TESTPEPT(Phospho (STY))IDE', 'TES(Phospho (STY))TPEPTIDE', 'T(Phospho (STY))ESTPEPTIDE']) == pd.Series(['TESTPEPTIDE.1.p1/p3/p8', 'TESTPEPTIDE'])).all()
    assert (transfer.get_modified_and_raw_sequence(['_S(Phospho (STY))RSAS(Phospho (STY))LRR_', '_SRS(Phospho (STY))AS(Phospho (STY))LRR_']) == pd.Series(['TESTPEPTIDE.1.p1/p3/p8', 'TESTPEPTIDE'])).all()


def test_generate_modified_sequence_annotation_ambiguous_raw_sequence():
    assert (transfer.get_modified_and_raw_sequence(['TESTPEPT(Phospho (STY))IDE', 'SET(Phospho (STY))TPEPTIDE', 'T(Phospho (STY))ESTPEPTIDE']) == pd.Series([np.nan, np.nan])).all()


def test_purge_mrc_files():
    assert False


def test_purge_modified_sequence():
    assert False


def test_clean_modified_sequence():
    assert False


def test_remove_nan_values():
    assert False


def test_get_single_item():
    assert False


def test_transfer():
    assert False


def test_generate_summary_file():
    assert False


def test_substitute_modifications():
    assert False


def test_flag_ambiguous_clusters():
    assert False


def test_get_main_object():
    assert False


def test_count_phos():
    assert False


def test_count_clustering_parameters():
    assert False


def test_cluster_mzml_files():
    assert False


# Creating dataframes from strings: https://towardsdatascience.com/67b0c2b71e6a
@pytest.fixture
def summary_missing_precursors():
    df_string = """new_type; Raw file; Charge;     m/z; Intensity; MS scan number
                 MULTI-MSMS;   file_1;      2;  392.12;   1238800;             12
                       MSMS;   file_1;      2; 1396.12;          ;             13
                       MSMS;   file_2;      1; 2396.12;          ;             14
                       MSMS;   file_2;      4; 3396.12;          ;             36
                       MSMS;   file_3;      3; 3396.12;          ;             36
                       MSMS;   file_4;      4; 3396.12;          ;             36"""
    return pd.read_csv(io.StringIO(df_string), delimiter=';', skipinitialspace=True)


@pytest.fixture
def allpeptides_missing_precursors():
    df_string = """new_type; Raw file; Charge;     m/z; Intensity; Min scan number; Max scan number
                      MULTI;   file_1;      2;  392.11;   1238800;              12;              15
                      MULTI;   file_1;      2; 1396.13;   2348800;              11;              19
                      MULTI;   file_2;      1; 2396.12;   3348800;               1;              12
                      MULTI;   file_3;      4; 3396.12;   4458800;              22;              45"""
    return pd.read_csv(io.StringIO(df_string), delimiter=';', skipinitialspace=True)

