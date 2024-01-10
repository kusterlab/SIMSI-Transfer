import pandas as pd
import pytest

import simsi_transfer.merging_functions as mf


class TestMergeWithMsmsTxt:
    @pytest.fixture
    def example_dataframes(self):
        # Create example dataframes for testing
        tempfile_data = {
            "Raw file": ["file1", "file2", "file3"],
            "scanID": [1, 2, 3],
            "temp_col": ["a", "b", "c"],
        }
        msms_data = {
            "Raw file": ["file1", "file2", "file4"],
            "scanID": [1, 2, 4],
            "msms_col": ["x", "y", "z"],
        }
        tempfile = pd.DataFrame(tempfile_data)
        msms_df = pd.DataFrame(msms_data)
        return tempfile, msms_df

    def test_merge_with_msmstxt(self, example_dataframes):
        tempfile, msms_df = example_dataframes
        merged_df = mf.merge_with_msmstxt(tempfile, msms_df)

        # Check if the original columns are present
        assert "temp_col" in merged_df.columns
        assert "msms_col" in merged_df.columns

        # Check specific values after merging
        assert merged_df.loc[0, "temp_col"] == "a"
        assert merged_df.loc[1, "msms_col"] == "y"

        # Check for NaN values where there is no match
        assert pd.isna(merged_df.loc[2, "msms_col"])

        # Check if the length is correct
        assert len(merged_df) == max(len(tempfile), len(msms_df))


class TestMergeSummaryWithEvidence:

    @pytest.fixture
    def example_dataframes(self):
        # Create example dataframes for testing
        summary_data = {'Modified sequence': ['AAA', 'BBB', 'CCC'],
                        'Raw file': ['file1', 'file2', 'file3'],
                        'Charge': [2, 3, 2],
                        'Other_column': ['x', 'y', 'z']}
        evidence_data = {'evidence_ID': [1, 2, 3],
                         'Modified sequence': ['AAA', 'BBB', 'DDD'],
                         'Raw file': ['file1', 'file2', 'file4'],
                         'Charge': [2, 3, 2],
                         'Intensity': [100, 200, 150],
                         'Leading proteins': ['protein1', 'protein2', 'protein3'],
                         'Type': ['type1', 'type2', 'type3'],
                         'Calibrated retention time': [10.5, 11.2, 9.8],
                         'Calibrated retention time start': [10.0, 11.0, 9.0],
                         'Calibrated retention time finish': [11.0, 12.0, 10.0],
                         'Retention time calibration': ['calib1', 'calib2', 'calib3']}
        summary = pd.DataFrame(summary_data)
        evidence = pd.DataFrame(evidence_data)
        return summary, evidence

    def test_merge_summary_with_evidence(self, example_dataframes):
        # Arrange
        summary, evidence = example_dataframes

        # Act
        merged_summary = mf.merge_summary_with_evidence(summary, evidence)

        # Assert
        # Check if the merged dataframe has the expected columns
        expected_columns = ['Modified sequence', 'Raw file', 'Charge', 'Other_column',
                             'evidence_ID', 'Intensity', 'Leading proteins', 'Type',
                             'Calibrated retention time', 'Calibrated retention time start',
                             'Calibrated retention time finish', 'Retention time calibration']
        assert all(column in merged_summary.columns for column in expected_columns)

        print(merged_summary)

        # Check if the merged dataframe has the correct values
        assert merged_summary.loc[0, 'Intensity'] == 100
        assert merged_summary.loc[1, 'evidence_ID'] == 2
        assert pd.isna(merged_summary.loc[2, 'Intensity'])  # Check NaN value for non-matching row

