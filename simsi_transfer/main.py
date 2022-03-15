import sys
import os
import warnings
from pathlib import Path
from datetime import datetime
import logging

import pandas as pd

from .IO_functions import export_summary_file, open_msms_txt, open_msmsscans_txt, open_maracluster_clusters, \
    parse_args, open_summary_txt, open_evidence_txt, export_simsi_evidence_file
from .processing_functions import generate_summary_file, flag_ambiguous_clusters, transfer, \
    count_clustering_parameters, count_phos, build_evidence, remove_unidentified_scans
from . import thermo_raw as raw
from . import maracluster as cluster
from . import tmt_extractor

logger = logging.getLogger(__name__)


def main(argv):
    mq_txt_folder, raw_folder, pvals, output_folder, num_threads, tmt_correction_file, ms_level = parse_args(argv)

    starttime = datetime.now()

    mq_txt_folder, raw_folder, output_folder, tmt_correction_file = Path(mq_txt_folder), Path(raw_folder), Path(
        output_folder), Path(tmt_correction_file)

    logger.info(f'Input parameters:')
    logger.info(f"MaxQuant txt folder = {mq_txt_folder}")
    logger.info(f"Raw file folder = {raw_folder}")
    logger.info(f"Stringencies = {','.join(map(str, pvals))}")
    logger.info(f"Output folder = {output_folder}")
    logger.info(f"Number of threads = {num_threads}")
    logger.info(f"TMT correction file = {tmt_correction_file}")
    logger.info(f"TMT MS level = {ms_level}")
    logger.info('')

    logger.info(f'Starting SIMSI-Transfer')
    logger.info('')

    if not output_folder.is_dir():
        output_folder.mkdir(parents=True)

    logger.info(f'Converting .raw files')
    mzml_folder = output_folder / Path('mzML')
    mzml_files = raw.convert_raw_mzml_batch(raw_folder, mzml_folder, num_threads)

    logger.info(f'Clustering .mzML files')
    maracluster_folder = output_folder / Path('maracluster_output')
    cluster.cluster_mzml_files(mzml_files, pvals, maracluster_folder, num_threads)

    logger.info(f'Reading in MaxQuant msmsscans.txt file')
    msmsscanstxt, tmt = open_msmsscans_txt(mq_txt_folder, ms_level)
    logger.info(msmsscanstxt.columns)

    if ms_level == 'ms3':
        logger.info(f'Extracting correct reporter ion intensities from .mzML files')
        extracted_folder = output_folder / Path('extracted')
        tmt_extractor.extract_tmt_reporters(mzml_files, extracted_folder, tmt_correction_file, num_threads)

        corrected_tmt = pd.DataFrame(
            columns=['raw_file', 'scanID', 'raw_TMT1', 'raw_TMT2', 'raw_TMT3', 'raw_TMT4', 'raw_TMT5', 'raw_TMT6',
                     'raw_TMT7', 'raw_TMT8', 'raw_TMT9', 'raw_TMT10', 'raw_TMT11', 'raw_TMT12', 'raw_TMT13',
                     'corr_TMT1',
                     'corr_TMT2', 'corr_TMT3', 'corr_TMT4', 'corr_TMT5', 'corr_TMT6', 'corr_TMT7', 'corr_TMT8',
                     'corr_TMT9', 'corr_TMT10', 'corr_TMT11'])

        for file in os.listdir(extracted_folder):
            file = pd.read_csv(Path(extracted_folder / file), sep='\t')
            corrected_tmt = corrected_tmt.append(file)
        corrected_tmt = corrected_tmt.reset_index(drop=True)
        corrected_tmt = corrected_tmt.rename(columns={
            'raw_file': 'Raw file',
            'raw_TMT1': 'Reporter intensity 1',
            'raw_TMT2': 'Reporter intensity 2',
            'raw_TMT3': 'Reporter intensity 3',
            'raw_TMT4': 'Reporter intensity 4',
            'raw_TMT5': 'Reporter intensity 5',
            'raw_TMT6': 'Reporter intensity 6',
            'raw_TMT7': 'Reporter intensity 7',
            'raw_TMT8': 'Reporter intensity 8',
            'raw_TMT9': 'Reporter intensity 9',
            'raw_TMT10': 'Reporter intensity 10',
            'raw_TMT11': 'Reporter intensity 11',
            'corr_TMT1': 'Reporter intensity corrected 1',
            'corr_TMT2': 'Reporter intensity corrected 2',
            'corr_TMT3': 'Reporter intensity corrected 3',
            'corr_TMT4': 'Reporter intensity corrected 4',
            'corr_TMT5': 'Reporter intensity corrected 5',
            'corr_TMT6': 'Reporter intensity corrected 6',
            'corr_TMT7': 'Reporter intensity corrected 7',
            'corr_TMT8': 'Reporter intensity corrected 8',
            'corr_TMT9': 'Reporter intensity corrected 9',
            'corr_TMT10': 'Reporter intensity corrected 10',
            'corr_TMT11': 'Reporter intensity corrected 11'})
        corrected_tmt['Raw file'] = corrected_tmt['Raw file'].str.replace(r'.mzML$', '')
        logger.info(corrected_tmt.sample(10))

        msmsscanstxt = pd.merge(left=msmsscanstxt, right=corrected_tmt, on=['Raw file', 'scanID'], how='left',
                                validate='one_to_one')
        logger.info(msmsscanstxt.sample(50))

    logger.info(f'Reading in MaxQuant msms.txt file and filtering out decoy hits')
    msmstxt = open_msms_txt(mq_txt_folder)
    # TODO: check if we should also transfer decoys
    msmstxt = msmstxt[msmstxt['Reverse'] != '+']

    logger.info(f'Reading in MaxQuant evidence.txt file and filtering out decoy hits')
    evidencetxt = open_evidence_txt(mq_txt_folder)
    # TODO: check if we should also transfer decoys
    evidencetxt = evidencetxt[evidencetxt['Reverse'] != '+']

    logger.info(f'Reading in MaxQuant summary.txt file')
    summarytxt = open_summary_txt(mq_txt_folder)

    statistics = dict()

    for pval in ['p' + str(i) for i in pvals]:
        logger.info('')

        logger.info(f'Starting MaxQuant and MaRaCluster file merge for {pval}.')

        clusterfile = open_maracluster_clusters(maracluster_folder, pval)
        summary = generate_summary_file(msmsscanstxt, msmstxt, summarytxt, clusterfile)
        del clusterfile

        export_summary_file(summary, output_folder, pval, state='merged')
        logger.info(f'Finished file merge.')

        logger.info(f'Starting cluster-based identity transfer for {pval}.')
        summary['phosphogroups'] = summary['Modified sequence'].apply(count_phos)
        summary = flag_ambiguous_clusters(summary)
        summary = transfer(summary)
        export_summary_file(summary, output_folder, pval, state='transferred')
        logger.info(f'Finished identity transfer.')

        logger.info(f'Building msms.csv file for {pval}.')
        summary = remove_unidentified_scans(summary)
        export_summary_file(summary, output_folder, pval, state='filtered')
        logger.info(f'Finished msms.csv assembly.')

        statistics[pval] = count_clustering_parameters(summary)

        logger.info(f'Starting evidence_transferred.txt building for {pval}.')
        evidence = build_evidence(summary, evidencetxt, tmt)
        export_simsi_evidence_file(evidence, output_folder, pval)
        logger.info(f'Finished evidence_transferred.txt building.')
        logger.info('')

    endtime = datetime.now()
    logger.info(f'Successfully finished transfers for all stringencies.')
    logger.info('')
    logger.info(f"SIMSI-Transfer finished in {(endtime - starttime).total_seconds()} seconds (wall clock).")


if __name__ == '__main__':
    main(sys.argv[1:])
