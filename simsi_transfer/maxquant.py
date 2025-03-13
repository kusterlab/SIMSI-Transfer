import logging
import re
from pathlib import Path

import pandas as pd
import numpy as np


logger = logging.getLogger(__name__)


def get_plex(input_folders):
    columns = pd.read_csv(
        input_folders[0] / Path("msms.txt"), nrows=0, sep="\t"
    ).columns.tolist()
    substring = re.compile(r"^Reporter intensity (\d{1,2})$")
    reporters = [i for i in columns if re.match(substring, i)]
    plex_number = max([int(re.search(substring, i).group(1)) for i in reporters])
    logger.info(f"Automatic plex detection: {plex_number}plex detected.")
    return plex_number


def read_msmsscans_txt(mq_txt_folder, tmt_requantify, plex):
    """
    Open msms.txt output file and subselect relevant columns
    :param mq_txt_folder: Processing path containing the 'combined' folder from MQ search
    :param tmt_requantify: Boolean indicating whether TMT reporter intensities need to be corrected later on
    :param plex: Number of TMT channels used in the experiment
    :return: truncated msmsscans.txt dataframe
    """
    columns = {
        "Raw file": "object",
        "Scan number": "int32",
        "m/z": "float32",
        "Mass": "float32",
        "Retention time": "float32",
        "Precursor full scan number": "int32",
        "MS scan number": "int32",
    }
    if not tmt_requantify:
        columns |= {f"Reporter intensity {i}": "float32" for i in range(1, plex + 1)}
        columns |= {
            f"Reporter intensity corrected {i}": "float32" for i in range(1, plex + 1)
        }
    msmsscans = pd.read_csv(
        mq_txt_folder / Path("msmsScans.txt"),
        sep="\t",
        usecols=columns.keys(),
        dtype=columns,
        engine="pyarrow",
    )

    msmsscans = msmsscans.rename(columns={"Scan number": "scanID"})
    return msmsscans


def read_msms_txt(mq_txt_folder):
    """
    Open msms.txt output file and subselect relevant columns
    :param mq_txt_folder: Processing path containing the 'combined' folder from MQ search
    :return: truncated msms.txt dataframe
    """
    columns = {
        "Raw file": "object",
        "Scan number": "int32",
        "Sequence": "object",
        "Modified sequence": "object",
        "Phospho (STY) Probabilities": "object",
        "Length": "int8",
        "Modifications": "object",
        "Missed cleavages": "int8",
        "Proteins": "object",
        "Gene Names": "object",
        "Protein Names": "object",
        "Charge": "int8",
        "Mass error [ppm]": "float32",
        "PIF": "float32",
        "Precursor Intensity": "float32",
        "PEP": "float32",
        "Score": "float32",
        "Delta score": "float32",
        "Reverse": "category",
    }

    columns_present = pd.read_csv(
        mq_txt_folder / Path("msms.txt"), nrows=0, sep="\t"
    ).columns.tolist()

    columns_to_read = {c: t for c, t in columns.items() if c in columns_present}

    msmstxt = pd.read_csv(
        mq_txt_folder / Path("msms.txt"),
        sep="\t",
        usecols=columns_to_read.keys(),
        dtype=columns_to_read,
        engine="pyarrow",
    )

    for col, dtype in columns.items():
        if col not in msmstxt.columns:
            logger.warning(f"Missing column in msms.txt, filled with numpy NaN: {col}")
            msmstxt[col] = np.NaN
            msmstxt[col] = msmstxt[col].astype(dtype)

    msmstxt = msmstxt.rename(columns={"Scan number": "scanID"})
    return msmstxt


def read_evidence_txt(mq_txt_folder):
    """
    Open msms.txt output file and subselect relevant columns
    :param mq_txt_folder: Processing path containing the 'combined' folder from MQ search
    :return: truncated evidence.txt dataframe
    """
    columns = {
        "Sequence": "object",
        "Modified sequence": "object",
        "Leading proteins": "object",
        "Raw file": "object",
        "Experiment": "object",
        "Fraction": "int8",
        "Charge": "int8",
        "Calibrated retention time": "float32",
        "Retention time": "float32",
        "Retention length": "float32",
        "Calibrated retention time start": "float32",
        "Calibrated retention time finish": "float32",
        "Retention time calibration": "float32",
        "Type": "object",
        "Intensity": "float32",
        "Reverse": "category",
    }

    columns_present = pd.read_csv(
        mq_txt_folder / Path("evidence.txt"), nrows=0, sep="\t"
    ).columns.tolist()

    columns_to_read = {c: t for c, t in columns.items() if c in columns_present}

    evidence = pd.read_csv(
        mq_txt_folder / Path("evidence.txt"),
        sep="\t",
        usecols=columns_to_read.keys(),
        dtype=columns_to_read,
        engine="pyarrow",
    )
    return evidence


def read_allpeptides_txt(mq_txt_folder):
    """
    Open msms.txt output file and subselect relevant columns
    :param mq_txt_folder: Processing path containing the 'combined' folder from MQ search
    :return: truncated evidence.txt dataframe
    """
    columns = {
        "Raw file": "object",
        "Type": "object",
        "Charge": "int8",
        "m/z": "float32",
        "Retention time": "float32",
        "Retention length": "float32",
        "Min scan number": "int32",
        "Max scan number": "int32",
        "Intensity": "float32",
    }
    allpeptides = pd.read_csv(
        mq_txt_folder / Path("allPeptides.txt"),
        sep="\t",
        usecols=columns.keys(),
        dtype=columns,
        engine="pyarrow",
    )
    return allpeptides


def fill_missing_min_max_scans(allpeptides, msms):
    if allpeptides[['Min scan number', 'Max scan number']].isna().any().any():
        msms_max_scans = msms.groupby('Raw file')['Precursor full scan number'].max()
        allpeptides['Max scan number'].fillna(allpeptides['Raw file'].map(msms_max_scans), inplace=True)
        allpeptides['Min scan number'].fillna(1, inplace=True)
        allpeptides['Max scan number'] = allpeptides['Max scan number'].astype(int)
        allpeptides['Min scan number'] = allpeptides['Min scan number'].astype(int)
    return allpeptides


def get_rawfile_metadata(evidence_txt):
    metadata_columns = ["Raw file", "Experiment", "Fraction"]
    metadata_columns_available = [
        x for x in metadata_columns if x in evidence_txt.columns
    ]
    metadata_columns_unavailable = [
        x for x in metadata_columns if x not in evidence_txt.columns
    ]

    meta_df = evidence_txt[metadata_columns_available].drop_duplicates()
    meta_df[metadata_columns_unavailable] = 1
    return meta_df


if __name__ == "__main__":
    raise NotImplementedError("Do not run this script.")
