import pytest
import pandas as pd

import simsi_transfer.maxquant as mq

pd.set_option('display.max_columns', None)


@pytest.fixture
def allpeptides_df():
    return pd.DataFrame({
        'Raw file': ['file1', 'file2', 'file3', 'file4', 'file5'],
        'Min scan number': [1, 2, None, None, 5],
        'Max scan number': [101, None, 103, None, 105]
    })


@pytest.fixture
def msms_df():
    return pd.DataFrame({
        'Raw file': ['file1', 'file2', 'file3', 'file4', 'file5'],
        'Precursor full scan number': [111, 112, 113, 114, 115]
    })


def test_fill_missing_min_max_scans(allpeptides_df, msms_df):
    result = mq.fill_missing_min_max_scans(allpeptides_df, msms_df)
    expected = pd.DataFrame({
        'Raw file': ['file1', 'file2', 'file3', 'file4', 'file5'],
        'Min scan number': [1, 2, 1, 1, 5],
        'Max scan number': [101, 112, 103, 114, 105]
    })
    pd.testing.assert_frame_equal(result, expected)
