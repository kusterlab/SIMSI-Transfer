import sys
import os
from pathlib import Path
from datetime import datetime
import logging
import time

from . import command_line_interface as cli
from . import maxquant as mq
from . import simsi_output
from . import thermo_raw as raw
from . import maracluster as cluster
from . import tmt_processing
from . import transfer
from . import evidence
from . import version

__version__ = version.get_version_from_pyproject()
__copyright__ = '''Copyright (c) 2021-2022 Firas Hamood & Matthew The. All rights reserved. Written by Firas Hamood (firas.hamood@tum.de) and Matthew The (matthew.the@tum.de) at the Chair of Proteomics and Bioanalytics at the Technical University of Munich.'''

logger = logging.getLogger(__name__)


def main(argv):
    mq_txt_folders, raw_folders, pvals, output_folder, num_threads, tmt_correction_files, ms_level, tmt_requantify, \
        filter_decoys, ambiguity_decision, meta_input_file = cli.parse_args(argv)

    if not output_folder.is_dir():
        output_folder.mkdir(parents=True)

    module_name = ".".join(__name__.split(".")[:-1])
    file_logger = logging.FileHandler(output_folder / Path('SIMSI.log'))
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    formatter.converter = time.gmtime
    file_logger.setFormatter(formatter)
    logging.getLogger(module_name).addHandler(file_logger)

    starttime = datetime.now()

    logger.info(f'SIMSI-Transfer version {__version__}')
    logger.info(f'{__copyright__}')
    logger.info(f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}')

    logger.info(f'Input parameters:')
    logger.info(f"MaxQuant txt folder = {mq_txt_folders}")
    logger.info(f"Raw file folder = {raw_folders}")
    logger.info(f"Stringencies = {','.join(map(str, pvals))}")
    logger.info(f"Output folder = {output_folder}")
    logger.info(f"Number of threads = {num_threads}")
    logger.info(f"TMT correction file = {tmt_correction_files}")
    logger.info(f"TMT MS level = {ms_level}")
    logger.info('')

    logger.info(f'Starting SIMSI-Transfer')
    logger.info('')

    logger.info(f'Converting .raw files')
    mzml_folder = output_folder / Path('mzML')
    mzml_files = raw.convert_raw_mzml_batch(raw_folders, mzml_folder, num_threads)

    logger.info(f'Clustering .mzML files')
    cluster_result_folder = output_folder / Path('maracluster_output')
    cluster.cluster_mzml_files(mzml_files, pvals, cluster_result_folder, num_threads)

    logger.info(f'Reading in MaxQuant msmsscans.txt file')
    msmsscans_mq = mq.process_and_concat(mq_txt_folders, mq.read_msmsscans_txt, tmt_requantify=tmt_requantify)

    if tmt_requantify:
        logger.info(f'Extracting correct reporter ion intensities from .mzML files')
        extracted_folder = output_folder / Path('extracted')
        # TODO: support multiple TMT correction files
        tmt_processing.extract_tmt_reporters(mzml_files, extracted_folder, tmt_correction_files[0], num_threads)
        
        corrected_tmt = tmt_processing.assemble_corrected_tmt_table(mzml_files, extracted_folder)
        msmsscans_mq = tmt_processing.merge_with_corrected_tmt(msmsscans_mq, corrected_tmt)

    logger.info(f'Reading in MaxQuant msms.txt file and filtering out decoy hits')
    msms_mq = mq.process_and_concat(mq_txt_folders, mq.read_msms_txt)
    # TODO: check if we should also transfer decoys
    if filter_decoys:
        msms_mq = msms_mq[msms_mq['Reverse'] != '+']

    logger.info(f'Reading in MaxQuant evidence.txt file and filtering out decoy hits')
    evidence_mq = mq.process_and_concat(mq_txt_folders, mq.read_evidence_txt)
    if filter_decoys:
        evidence_mq = evidence_mq[evidence_mq['Reverse'] != '+']
    rawfile_metadata = mq.get_rawfile_metadata(evidence_mq)

    logger.info(f'Reading in MaxQuant allPeptides.txt file')
    allpeptides_mq = mq.process_and_concat(mq_txt_folders, mq.read_allpeptides_txt)

    statistics = dict()

    for pval in ['p' + str(i) for i in pvals]:
        logger.info('')

        logger.info(f'Starting MaxQuant and MaRaCluster file merge for {pval}.')

        cluster_results = cluster.read_cluster_results(cluster_result_folder, pval)
        annotated_clusters = simsi_output.annotate_clusters(msmsscans_mq, msms_mq, rawfile_metadata, cluster_results)
        del cluster_results

        simsi_output.export_annotated_clusters(annotated_clusters, output_folder, pval)
        logger.info(f'Finished file merge.')

        logger.info(f'Starting cluster-based identity transfer for {pval}.')
        annotated_clusters = transfer.flag_ambiguous_clusters(annotated_clusters)
        msmsscans_simsi = transfer.transfer(annotated_clusters, ambiguity_decision=ambiguity_decision)
        simsi_output.export_msmsscans(msmsscans_simsi, output_folder, pval)
        logger.info(f'Finished identity transfer.')

        logger.info(f'Building SIMSI-Transfer msms.txt file for {pval}.')
        msms_simsi = simsi_output.remove_unidentified_scans(msmsscans_simsi)
        simsi_output.export_msms(msms_simsi, output_folder, pval)
        logger.info(f'Finished SIMSI-Transfer msms.txt assembly.')

        statistics[pval] = simsi_output.count_clustering_parameters(msms_simsi)

        logger.info(f'Starting SIMSI-Transfer evidence.txt building for {pval}.')
        evidence_simsi = evidence.build_evidence(msms_simsi, evidence_mq, allpeptides_mq)
        simsi_output.export_simsi_evidence_file(evidence_simsi, output_folder, pval)
        logger.info(f'Finished SIMSI-Transfer evidence.txt building.')
        logger.info('')

    endtime = datetime.now()
    logger.info(f'Successfully finished transfers for all stringencies.')
    logger.info('')
    logger.info(f"SIMSI-Transfer finished in {(endtime - starttime).total_seconds()} seconds (wall clock).")


if __name__ == '__main__':
    main(sys.argv[1:])
