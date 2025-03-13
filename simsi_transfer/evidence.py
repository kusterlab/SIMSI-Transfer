import logging
from typing import List

import numpy as np
import pandas as pd

from .merging_functions import merge_summary_with_evidence
from .utils import utils

logger = logging.getLogger(__name__)


def assign_evidence_type(summary: pd.DataFrame, type_column_name: str = "new_type"):
    """
    assign the updated Type column by checking if retention time of MS2 scan is within its evidence precursor rt window
    """
    summary[type_column_name] = np.NaN
    msms_in_rt_window = (
        summary["Retention time"]
        >= summary["Calibrated retention time start"]
        - summary["Retention time calibration"]
    ) & (
        summary["Retention time"]
        <= summary["Calibrated retention time finish"]
        - summary["Retention time calibration"]
    )
    summary.loc[msms_in_rt_window, type_column_name] = "MULTI-MSMS"
    summary.loc[~msms_in_rt_window, type_column_name] = "MSMS"

    # clear out the retention time columns for MS2 scans without precursor
    del_columns = [
        "Type",
        "Calibrated retention time",
        "Calibrated retention time start",
        "Calibrated retention time finish",
    ]
    summary.loc[~msms_in_rt_window, del_columns] = np.NaN
    return summary


def assign_missing_precursors(summary: pd.DataFrame, allpeptides: pd.DataFrame):
    """
    find precursors in allPeptides.txt for transferred MS2 scans in runs where the
    (peptide, charge) combination was previously not identified currently very slow
    (~500 MS2 scans / second).
    """
    group_key = ["Raw file", "Charge"]
    allpeptides_grouped = {
        group: df_grouped for group, df_grouped in allpeptides.groupby(group_key)
    }

    missing_precursor = summary["new_type"] == "MSMS"
    summary.loc[missing_precursor, ["new_type", "Intensity"]] = (
        summary.loc[missing_precursor]
        .apply(
            lambda x: match_precursor_grouped(x, allpeptides_grouped, group_key),
            axis=1,
            result_type="expand",
        )
        .to_numpy()
    )  # need to convert to numpy array, see: https://stackoverflow.com/questions/69954697/

    num_new_assigned_precursor = (
        summary.loc[missing_precursor]["new_type"] == "MULTI-MSMS"
    ).sum()
    logger.debug(
        f"Assigned precursor to {num_new_assigned_precursor} out of {missing_precursor.sum()} unassigned MS2 spectra"
    )
    return summary


def get_ppm_diff(mz1, mz2):
    return np.abs(mz1 - mz2) / mz1 * 1e6


def match_precursor_grouped(
    msms_scan: pd.Series,
    allpeptides_grouped: pd.core.groupby.DataFrameGroupBy,
    group_key: List[str],
):
    if tuple(msms_scan[group_key]) not in allpeptides_grouped:
        return ["MSMS", np.NaN]
    allpeptides = allpeptides_grouped[tuple(msms_scan[group_key])]
    return match_precursor(msms_scan, allpeptides)


def match_precursor(
    msms_scan: pd.Series, allpeptides: pd.DataFrame, ppm_tol: float = 20.0
):
    precursors = allpeptides[
        (get_ppm_diff(allpeptides["m/z"], msms_scan["m/z"]) < ppm_tol)
        & (allpeptides["Min scan number"] <= msms_scan["MS scan number"])
        & (allpeptides["Max scan number"] >= msms_scan["MS scan number"])
    ]
    if len(precursors.index) == 0:
        return ["MSMS", np.NaN]
    return ["MULTI-MSMS", precursors["Intensity"].values[0]]


def remove_duplicate_msms(summary: pd.DataFrame):
    """
    remove duplicate MS2 scans due to (peptide, raw_file, charge) combinations matching multiple precursors
    prioritize the entry that matched a precursor (MULTI-MSMS) or has the lower retention time
    """
    summary = summary.sort_values(
        by=["new_type", "Calibrated retention time"], ascending=[False, True]
    )
    summary = summary.drop_duplicates(subset=["summary_ID"], keep="first")
    return summary


def fill_missing_evidence_ids(summary):
    """
    fill the evidence_ID column of MS2 scans without a precursor
    """
    missing_evidence_id = summary["evidence_ID"].isna()

    start_id = int(summary["evidence_ID"].max() + 1)
    end_id = start_id + missing_evidence_id.sum()

    summary.loc[missing_evidence_id, "evidence_ID"] = range(start_id, end_id)
    summary["evidence_ID"] = summary["evidence_ID"].astype(int)
    return summary


def assign_evidence_feature(
    summary: pd.DataFrame, evidence: pd.DataFrame, allpeptides: pd.DataFrame
):
    # store number of rows of summary dataframe to check if we have the same number after merging
    summary_length_before_processing = len(summary.index)

    summary = merge_summary_with_evidence(summary, evidence)
    summary = assign_evidence_type(summary)
    summary = remove_duplicate_msms(summary)
    summary = assign_missing_precursors(summary, allpeptides)
    summary = fill_missing_evidence_ids(summary)

    if not len(summary) == summary_length_before_processing:
        raise ValueError(
            f"Number of summary entries changed during evidence feature assembly!"
            f"\n{summary_length_before_processing} before assembly,\n{len(summary)} after assembly."
        )

    if not summary["summary_ID"].nunique() == len(summary):
        raise ValueError("Some summary_IDs were detected multiple times!")

    # TODO: check if all MULTI-MSMS entries are allocated to the correct PSMs

    summary = summary.sort_values(by=["Sequence", "Modified sequence"])
    summary = summary.drop(columns=["Type"])
    return summary


def calculate_evidence_columns(summary: pd.DataFrame, plex: int):
    # replacing zeros with NaNs to count later
    logger.debug("Assigned evidence features; calculating column values...")
    reps = [f"Reporter intensity {i}" for i in range(1, plex + 1)]
    reps.extend([f"Reporter intensity corrected {i}" for i in range(1, plex + 1)])
    summary[reps] = summary[reps].replace({0: np.nan})
    summary = summary.sort_values(
        by=["Sequence", "Modified sequence", "Raw file", "Charge"]
    ).reset_index(drop=True)

    # add semicolon to columns which will be concatenated. this allows us to use the
    # fast "sum" aggfunc instead of a slow custom string function
    concat_cols = ["Proteins", "Leading proteins", "scanID", "summary_ID"]
    summary[concat_cols] = summary[concat_cols].apply(
        lambda column: column.astype(str) + ";"
    )

    summary_grouped = summary.groupby("evidence_ID")

    # Column generation
    # TODO: Sort for keeping best hit at top; best hit definition needed
    # TODO: Check every column; what is needed and is something missing?

    evidence = summary_grouped.agg(
        **{
            "Sequence": pd.NamedAgg(
                column="Sequence", aggfunc="first"
            ),  # from msms.txt
            "Length": pd.NamedAgg(column="Length", aggfunc="first"),  # from msms.txt
            "Modifications": pd.NamedAgg(
                column="Modifications", aggfunc="first"
            ),  # from msms.txt
            "Modified sequence": pd.NamedAgg(
                column="Modified sequence", aggfunc="first"
            ),  # from msms.txt
            "Missed cleavages": pd.NamedAgg(
                column="Missed cleavages", aggfunc="first"
            ),  # from msms.txt
            "Proteins": pd.NamedAgg(column="Proteins", aggfunc="sum"),  # from msms.txt
            "Leading proteins": pd.NamedAgg(column="Leading proteins", aggfunc="sum"),
            # from evidence.txt, NaN if scan not matched to precursor in evidence.txt
            "Gene Names": pd.NamedAgg(
                column="Gene Names", aggfunc="first"
            ),  # from msms.txt
            "Protein Names": pd.NamedAgg(
                column="Protein Names", aggfunc="first"
            ),  # from msms.txt
            "Type": pd.NamedAgg(
                column="new_type", aggfunc="first"
            ),  # calculated by SIMSI-Transfer
            "Raw file": pd.NamedAgg(
                column="Raw file", aggfunc="first"
            ),  # from msms.txt
            "Fraction": pd.NamedAgg(column="Fraction", aggfunc="first"),
            # from evidence.txt using rawfile_metadata dictionary
            "Experiment": pd.NamedAgg(column="Experiment", aggfunc="first"),
            # from evidence.txt using rawfile_metadata dictionary
            "Charge": pd.NamedAgg(column="Charge", aggfunc="first"),  # from msms.txt
            "m/z": pd.NamedAgg(column="m/z", aggfunc="first"),  # from msmsScans.txt
            "Mass": pd.NamedAgg(column="Mass", aggfunc="first"),  # from msmsScans.txt
            "Mass error [ppm]": pd.NamedAgg(
                column="Mass error [ppm]", aggfunc="min"
            ),  # from msms.txt, NaN if Type=MSMS
            "Retention time": pd.NamedAgg(
                column="Retention time", aggfunc="first"
            ),  # from msmsScans.txt
            "PEP": pd.NamedAgg(column="PEP", aggfunc="min"),  # from msms.txt
            "MS/MS count": pd.NamedAgg(
                column="Sequence", aggfunc="size"
            ),  # calculated by SIMSI-Transfer
            "MS/MS all scan numbers": pd.NamedAgg(
                column="scanID", aggfunc="sum"
            ),  # calculated by SIMSI-Transfer
            "MS/MS scan number": pd.NamedAgg(
                column="scanID", aggfunc="first"
            ),  # calculated by SIMSI-Transfer
            "Score": pd.NamedAgg(column="Score", aggfunc="max"),  # from msms.txt
            "Delta score": pd.NamedAgg(
                column="Delta score", aggfunc="max"
            ),  # from msms.txt
            "Intensity": pd.NamedAgg(column="Intensity", aggfunc="sum"),
            # from evidence.txt, supplemented from allPeptides.txt
            **{
                f"Reporter intensity corrected {i}": pd.NamedAgg(
                    column=f"Reporter intensity corrected {i}", aggfunc="sum"
                )
                for i in range(1, plex + 1)
            },
            # from msmsScans.txt or parsed from mzML files by SIMSI-Transfer
            **{
                f"Reporter intensity {i}": pd.NamedAgg(
                    column=f"Reporter intensity {i}", aggfunc="sum"
                )
                for i in range(1, plex + 1)
            },  # from msmsScans.txt or parsed from mzML files by SIMSI-Transfer
            **{
                f"Reporter intensity count {i}": pd.NamedAgg(
                    column=f"Reporter intensity {i}", aggfunc="count"
                )
                for i in range(1, plex + 1)
            },  # calculated by SIMSI-Transfer
            "Reverse": pd.NamedAgg(column="Reverse", aggfunc="first"),  # from msms.txt
            "summary_ID": pd.NamedAgg(
                column="summary_ID", aggfunc="sum"
            ),  # assigned by SIMSI-Transfer
            "Transferred spectra count": pd.NamedAgg(
                column="identification", aggfunc="sum"
            ),  # calculated by SIMSI-Transfer
        }
    )

    for col in [
        "Proteins",
        "Leading proteins",
        "MS/MS all scan numbers",
        "MS/MS scan number",
        "summary_ID",
    ]:
        evidence[col] = evidence[col].apply(utils.csv_unique)

    evidence["Transferred spectra count"] = evidence[
        "Transferred spectra count"
    ].str.count("t")

    # evidence["Leading proteins"] = evidence["Leading proteins"].apply(utils.csv_list_unique)
    evidence["id"] = evidence.index  # evidence_ID

    # Checking for entries without "Fraction"; these are caused when one or more raw files in a maxquant search generate
    # no results in evidence.txt and therefore the metadata cannot be fetched
    no_fraction_rows = evidence["Fraction"].isna()
    no_metadata_rawfiles = evidence.loc[no_fraction_rows, "Raw file"].unique()
    if no_metadata_rawfiles.size > 0:
        logger.warning(
            f'No metadata fetched from evidence for {no_metadata_rawfiles}; autofilling with Fraction = -1 and Experiment = "UNKNOWN"'
        )
    evidence.loc[no_fraction_rows, "Experiment"] = "UNKNOWN"
    evidence.loc[no_fraction_rows, "Fraction"] = -1
    evidence = evidence.astype(
        {
            "Length": "int8",
            "Missed cleavages": "int8",
            "Fraction": "int8",
            "Charge": "int8",
            "MS/MS scan number": "int32",
        }
    )

    mask = evidence["Proteins"].str.contains("CON__").fillna(False)
    evidence.loc[mask, "Potential contaminant"] = "+"
    evidence.loc[~mask, "Potential contaminant"] = ""

    evidence = evidence.sort_values(
        by=["Sequence", "Modified sequence", "Raw file", "Charge"]
    )
    return evidence


def build_evidence_grouped(
    summary: pd.DataFrame,
    evidence: pd.DataFrame,
    allpeptides: pd.DataFrame,
    plex: int,
    num_threads: int,
):
    """
    Optimized merging function to group dataframes by 'Raw file', perform the merge per group,
    and concatenate the results into one dataframe with minimal memory usage.
    :param summary: intermediate summary dataframe derived by transferring
    :param evidence: evidence dataframe to build new evidence on
    :param allpeptides: allpeptides.txt dataframe to get all precursor information
    :param plex: number of TMT channels
    :param num_threads: number of threads to use
    """
    # logger.debug(f"summary memory usage: {utils.get_dataframe_size(summary)}")
    # logger.debug(f"evidence memory usage: {utils.get_dataframe_size(evidence)}")
    # logger.debug(f"allPeptides memory usage: {utils.get_dataframe_size(allpeptides)}")

    # Group dataframes by 'Raw file'
    summary_groups = summary.groupby("Raw file")
    evidence_groups = evidence.groupby("Raw file")
    allpeptides_groups = allpeptides.groupby("Raw file")

    multithreading = num_threads > 1
    if multithreading:
        from job_pool import JobPool

        job_pool = JobPool(
            processes=num_threads,
            write_progress_to_logger=True,
            total_jobs=len(summary_groups),
        )

    # Iterate through each group and merge
    evidences = []
    for raw_file, summary_group in summary_groups:
        if raw_file not in evidence_groups.groups:
            logger.warning(f"{raw_file} missing in evidence.txt, skipping this file")
            continue
        if raw_file not in allpeptides_groups.groups:
            logger.warning(f"{raw_file} missing in allPeptides.txt, skipping this file")
            continue

        evidence_group = evidence_groups.get_group(raw_file)
        allpeptides_group = allpeptides_groups.get_group(raw_file)

        args = (summary_group, evidence_group, allpeptides_group, plex)
        if multithreading:
            job_pool.applyAsync(build_evidence, args)
        else:
            evidences.append(build_evidence(*args))

    if multithreading:
        evidences = job_pool.checkPool()

    return pd.concat(evidences, ignore_index=True)


def build_evidence(
    summary: pd.DataFrame, evidence: pd.DataFrame, allpeptides: pd.DataFrame, plex: int
):
    evidence = evidence[evidence["Type"] != "MSMS"]
    evidence = evidence.sort_values(
        by=[
            "Sequence",
            "Modified sequence",
            "Raw file",
            "Calibrated retention time start",
        ]
    )
    evidence.insert(len(evidence.columns), "evidence_ID", range(len(evidence)))
    summary = assign_evidence_feature(summary, evidence, allpeptides)
    evidence = calculate_evidence_columns(summary, plex)
    return evidence
