import pandas as pd
import pytest
import numpy as np
import simsi_transfer.evidence as ev


class TestCalculateEvidenceColumns:
    @pytest.fixture
    def summary(self):
        # Create example dataframes for testing
        summary_data = {
            "evidence_ID": [1, 2, 1],
            "Sequence": ["AAA", "BBB", "AAA"],
            "Modifications": ["mod1", "mod2", "mod1"],
            "Modified sequence": ["MOD_AAA", "MOD_BBB", "MOD_AAA"],
            "Raw file": ["file1", "file2", "file1"],
            "Fraction": [1, 2, 1],
            "Experiment": ["Exp1", "Exp2", "Exp1"],
            "m/z": [123.45, 234.56, 123.45],
            "Mass": [1000, 1500, 1000],
            "Mass error [ppm]": [1.5, 2.0, np.nan],
            "Retention time": [50.0, 60.0, 70.0],
            "Charge": [2, 3, 2],
            "Length": [10, 15, 12],
            "Missed cleavages": [1, 2, 0],
            "Proteins": ["ProteinA;ProteinC", "ProteinB", "ProteinC"],
            "Leading proteins": ["LeadA", "LeadB", "LeadC;LeadA"],
            "Gene Names": ["GeneA", "GeneB", "GeneC"],
            "Protein Names": ["NameA", "NameB", "NameC"],
            "new_type": ["TypeA", "TypeB", "TypeC"],
            "Intensity": [100, 150, 200],
            "Reporter intensity 1": [10, 15, 20],
            "Reporter intensity 2": [5, 7, 10],
            "Reporter intensity corrected 1": [8, 12, 18],
            "Reporter intensity corrected 2": [4, 6, 9],
            "Reverse": ["", "+", ""],
            "scanID": [101, 102, 103],
            "PEP": [0.01, 0.02, 0.03],
            "Score": [50, 60, 70],
            "Delta score": [5, 6, 7],
            "summary_ID": [1, 2, 3],
            "identification": ["t", "t", "d"],
        }
        return pd.DataFrame(summary_data)

    def test_calculate_evidence_columns(self, summary):
        plex = 2  # Assuming a 2-plex experiment, adjust as needed

        # Act
        calculated_evidence = ev.calculate_evidence_columns(summary, plex)

        # Assert
        # Check if the calculated evidence dataframe has the expected columns
        expected_columns = [
            "id",
            "Sequence",
            "Length",
            "Modifications",
            "Modified sequence",
            "Missed cleavages",
            "Proteins",
            "Leading proteins",
            "Gene Names",
            "Protein Names",
            "Type",
            "Raw file",
            "Fraction",
            "Experiment",
            "Charge",
            "m/z",
            "Mass",
            "Mass error [ppm]",
            "Retention time",
            "PEP",
            "MS/MS count",
            "MS/MS all scan numbers",
            "MS/MS scan number",
            "Score",
            "Delta score",
            "Intensity",
            "Reporter intensity corrected 1",
            "Reporter intensity corrected 2",
            "Reporter intensity 1",
            "Reporter intensity 2",
            "Reporter intensity count 1",
            "Reporter intensity count 2",
            "Reverse",
            "summary_ID",
            "Transferred spectra count",
        ]
        assert all(column in calculated_evidence.columns for column in expected_columns)

        # Check if the calculated evidence dataframe has the correct values
        assert calculated_evidence.loc[1, "Intensity"] == 300
        assert calculated_evidence.loc[1, "Fraction"] == 1
        assert calculated_evidence.loc[1, "Proteins"] == "ProteinA;ProteinC"
        assert calculated_evidence.loc[1, "Leading proteins"] == "LeadA;LeadC"
        assert calculated_evidence.loc[1, "Transferred spectra count"] == 1
        assert calculated_evidence.loc[1, "Reporter intensity corrected 1"] == 26
        assert calculated_evidence.loc[1, "Reporter intensity count 1"] == 2
        assert calculated_evidence.loc[1, "summary_ID"] == "1;3"
        
        assert calculated_evidence.loc[2, "Reverse"] == "+"
        assert calculated_evidence.loc[2, "Transferred spectra count"] == 1
        
