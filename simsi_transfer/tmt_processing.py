import sys
import re
from typing import List
from pathlib import Path
import logging

import pandas as pd
import numpy as np
from pyteomics import mzml

from .utils import utils

logger = logging.getLogger(__name__)


def get_tmt_columns(plex):
    return [f"raw_TMT{i}" for i in range(1, plex + 3)], [
        f"corr_TMT{i}" for i in range(1, plex + 1)
    ]


def get_correction_factors(correction_factor_path: Path, plex_size: int):
    # correction = np.array([[100, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 126 C Tag
    #                        [0.0, 100, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 127 N Tag
    #                        [0.0, 0.0, 100, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 127 C Tag
    #                        [0.0, 0.0, 0.0, 100, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 128 N Tag
    #                        [0.0, 0.0, 0.0, 0.0, 100, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 128 C Tag
    #                        [0.0, 0.0, 0.0, 0.0, 0.0, 100, 0.0, 0.0, 0.0, 0.0, 0.0],  # 129 N Tag
    #                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100, 0.0, 0.0, 0.0, 0.0],  # 129 C Tag
    #                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100, 0.0, 0.0, 0.0],  # 130 N Tag
    #                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100, 0.0, 0.0],  # 130 C Tag
    #                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100, 0.0],  # 131 N Tag
    #                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100],  # 131 C Tag
    #                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 132 N Overflow
    #                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # 132 C Overflow
    #                        ])
    correction = np.zeros(shape=(plex_size + 2, plex_size))
    for i in range(correction.shape[1]):
        correction[i, i] = 100

    # Theoretical TMT Masses in m/z; same for standard TMT and TMTpro, different for sixplex though!
    all_tmt_masses = np.array(
        [
            126.127726,
            127.124761,
            127.131081,
            128.128116,
            128.134436,
            129.131471,
            129.137790,
            130.134825,
            130.141145,
            131.138180,
            131.144499,
            132.141535,
            132.147855,
            133.144890,
            133.151210,
            134.148245,
            134.154565,
            135.151600,
        ]
    )

    tmt_masses = all_tmt_masses[: plex_size + 2]
    if plex_size == 6:
        tmt_masses = np.array([all_tmt_masses[i] for i in [0, 1, 4, 5, 8, 9, 10, 11]])

    np.set_printoptions(linewidth=200)
    if correction_factor_path.is_file():
        correction_dataframe = pd.read_csv(correction_factor_path, sep="\t")

        for i in range(plex_size):
            if i not in [0, 1, 2, 3]:
                correction[i - 4, i] = correction_dataframe.iloc[i][
                    "Correction factor -2 [%]"
                ]
            if i not in [0, 1]:
                correction[i - 2, i] = correction_dataframe.iloc[i][
                    "Correction factor -1 [%]"
                ]
            correction[i + 2, i] = correction_dataframe.iloc[i][
                "Correction factor +1 [%]"
            ]
            if i not in [plex_size - 1, plex_size - 2, plex_size - 3]:
                correction[i + 4, i] = correction_dataframe.iloc[i][
                    "Correction factor +2 [%]"
                ]
    elif str(correction_factor_path) != ".":
        logger.warning(
            f"Could not find reporter ion correction file at {str(correction_factor_path)}, no correction factors will be applied"
        )

    # Normalize correction factors
    correction_normalized = correction / correction.sum(axis=0)
    return tmt_masses, correction_normalized


def extract_tmt_reporters(
    mzml_files: List[Path],
    output_path: Path,
    correction_factor_paths: List[Path],
    plex: int,
    extraction_level: int,
    num_threads: int = 1,
):
    """
    Takes about 1.5 minute for a 700MB file with 40k MS2 scans
    """
    if not output_path.is_dir():
        output_path.mkdir(parents=True)

    if num_threads > 1:
        from job_pool import JobPool

        processing_pool = JobPool(processes=num_threads, write_progress_to_logger=True)

    for mzml_file, correction_factor_path in zip(mzml_files, correction_factor_paths):
        args = (mzml_file, output_path, correction_factor_path, extraction_level, plex)
        if num_threads > 1:
            processing_pool.applyAsync(extract_and_correct_reporters, args)
        else:
            extract_and_correct_reporters(*args)

    if num_threads > 1:
        processing_pool.checkPool()


def extract_and_correct_reporters(
    mzml_file: Path,
    output_path: str,
    correction_factor_path: Path,
    extraction_level: int,
    plex: int,
):
    tmt_masses, correction_normalized = get_correction_factors(
        correction_factor_path, plex_size=plex
    )

    tolerance = 6 * 1e-3 / 2
    tmt_upper = tmt_masses + tolerance
    tmt_lower = tmt_masses - tolerance

    tmt_raw_columns, tmt_corrected_columns = get_tmt_columns(plex)

    convert_dict_raw = {k: "float32" for k in tmt_raw_columns}
    convert_dict_corr = {k: "float32" for k in tmt_corrected_columns}
    convert_dict_other = {"raw_file": "str", "scanID": "int32"}
    convert_dict = {**convert_dict_raw, **convert_dict_corr, **convert_dict_other}
    dfcol = convert_dict.keys()

    output_file = get_extracted_tmt_file_name(output_path, mzml_file)
    if Path(output_file).is_file():
        logger.debug(
            f"Found extracted reporter ions at {output_file}, skipping extraction"
        )
        return

    logger.info("Performing extraction for " + mzml_file.name)
    fileframe = pd.DataFrame(columns=dfcol)

    seen_scan_ids = set()
    with mzml.read(str(mzml_file)) as reader:
        for i, item in enumerate(reader):
            if item["ms level"] != extraction_level:
                continue

            scanseries = pd.Series(index=dfcol, dtype="float32")

            if extraction_level == 2:
                scanseries["scanID"] = re.search(r"scan=(\d+)", item["id"])[1]
            else:
                # supposed to find parent MS2 spectrum for MS3 by looking into precursorList/precursor/spectrumRef
                scanseries["scanID"] = re.search(
                    r"scan=(\d+)", item["precursorList"]["precursor"][0]["spectrumRef"]
                )[1]

            if scanseries["scanID"] in seen_scan_ids:
                logger.warning(
                    f"Found duplicate MS3 spectrum for MS2 spectrum with scan number {scanseries['scanID']} in {mzml_file}, known bug in ThermoRawFileParser..."
                )
                continue
            seen_scan_ids.add(scanseries["scanID"])

            mz = np.array(item["m/z array"])
            intensity = np.array(item["intensity array"])
            for c, (low, upp) in enumerate(zip(tmt_lower, tmt_upper)):
                start_idx = int(np.searchsorted(mz, low))
                end_idx = int(np.searchsorted(mz, upp))
                scanseries[f"raw_TMT{c + 1}"] = intensity[start_idx:end_idx].sum()
            fileframe = pd.concat(
                [fileframe, scanseries.to_frame().T], ignore_index=True
            )
    fileframe["raw_file"] = mzml_file.name

    fileframe = fileframe.astype(convert_dict)

    # TMT correction
    logger.info("Extraction done, correcting TMT reporters for " + mzml_file.name)
    fileframe[tmt_corrected_columns] = pd.DataFrame(
        fileframe[tmt_raw_columns]
        .apply(
            lambda tmt: np.linalg.lstsq(correction_normalized, tmt, rcond=None)[
                0
            ].round(2),
            axis=1,
        )
        .tolist(),
        columns=tmt_corrected_columns,
        index=fileframe[tmt_corrected_columns].index,
    )
    # PANDAS WHERE! This retains the checked value if the condition is met, and replaces it where it is not met!
    fileframe[tmt_corrected_columns] = fileframe[tmt_corrected_columns].where(
        fileframe[tmt_corrected_columns] > 10, 0
    )
    fileframe.to_csv(output_file, sep="\t", index=False)


def get_extracted_tmt_file_name(output_path: str, mzml_file: Path):
    return f"{output_path}/ext_{mzml_file.name}.txt"


def read_extracted_tmt_file(extracted_tmt_file: Path, plex: int):
    columns = {
        "raw_file": "object",
        "scanID": "int32",
        **{f"raw_TMT{i}": "float32" for i in range(1, plex + 1)},
        **{f"corr_TMT{i}": "float32" for i in range(1, plex + 1)},
    }
    # engine="pyarrow" is not faster for these small files.
    return pd.read_csv(
        extracted_tmt_file,
        sep="\t",
        usecols=columns.keys(),
        dtype=columns,
    )


def assemble_corrected_tmt_table(
    mzml_files: List[Path], extracted_folder: Path, plex: int
):
    logger.info("Assembling corrected reporter ion tables")

    extracted_tmt_files = [
        get_extracted_tmt_file_name(extracted_folder, mzml_file)
        for mzml_file in mzml_files
    ]
    corrected_tmt = utils.process_and_concat(
        extracted_tmt_files, read_extracted_tmt_file, plex=plex
    )

    corrected_tmt = corrected_tmt.reset_index(drop=True)
    corrected_tmt = corrected_tmt.rename(
        columns={
            "raw_file": "Raw file",
            **{f"raw_TMT{i}": f"Reporter intensity {i}" for i in range(1, plex + 1)},
            **{
                f"corr_TMT{i}": f"Reporter intensity corrected {i}"
                for i in range(1, plex + 1)
            },
        }
    )
    corrected_tmt["Raw file"] = corrected_tmt["Raw file"].str.replace(
        r".mzML$", "", regex=True
    )
    return corrected_tmt


def merge_with_corrected_tmt(msmsscans_df: pd.DataFrame, corrected_tmt: pd.DataFrame):
    logger.info("Merging corrected reporter ion tables into msmsScans.txt")
    return pd.merge(
        left=msmsscans_df,
        right=corrected_tmt,
        on=["Raw file", "scanID"],
        how="left",
        validate="one_to_one",
    )


if __name__ == "__main__":
    input_files_arg = Path(sys.argv[1])
    extraction_level_arg = sys.argv[2]
    output_path_arg = Path(sys.argv[3])

    extract_tmt_reporters(
        [input_files_arg],
        output_path_arg,
        extraction_level=extraction_level_arg,
        plex=11,
    )
