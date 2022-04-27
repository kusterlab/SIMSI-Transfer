import sys
from pathlib import Path
from datetime import datetime
import logging

from . import command_line_interface as cli
from . import maxquant as mq
from . import simsi_output
from . import thermo_raw as raw
from . import maracluster as cluster
from . import tmt_processing
from . import transfer
from . import evidence

logger = logging.getLogger(__name__)


def main(argv):
    mq_txt_folder, raw_folder, pvals, output_folder, num_threads, tmt_correction_file, ms_level, tmt_requantify = cli.parse_args(
        argv)

    if not output_folder.is_dir():
        output_folder.mkdir(parents=True)

    module_name = ".".join(__name__.split(".")[:-1])
    logging.getLogger(module_name).addHandler(
        logging.FileHandler(output_folder / Path('SIMSI.log')))
    starttime = datetime.now()

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

    logger.info(f'Converting .raw files')
    mzml_folder = output_folder / Path('mzML')
    mzml_files = raw.convert_raw_mzml_batch(raw_folder, mzml_folder, num_threads)

    logger.info(f'Clustering .mzML files')
    cluster_result_folder = output_folder / Path('maracluster_output')
    cluster.cluster_mzml_files(mzml_files, pvals, cluster_result_folder, num_threads)

    logger.info(f'Reading in MaxQuant msmsscans.txt file')
    msmsscans_mq, tmt = mq.read_msmsscans_txt(mq_txt_folder, tmt_requantify)

    if tmt_requantify:
        logger.info(f'Extracting correct reporter ion intensities from .mzML files')
        extracted_folder = output_folder / Path('extracted')
        tmt_processing.extract_tmt_reporters(mzml_files, extracted_folder, tmt_correction_file, num_threads)
        corrected_tmt = tmt_processing.assemble_corrected_tmt_table(extracted_folder)

        msmsscans_mq = tmt_processing.merge_with_corrected_tmt(msmsscans_mq, corrected_tmt)

    logger.info(f'Reading in MaxQuant msms.txt file and filtering out decoy hits')
    msms_mq = mq.read_msms_txt(mq_txt_folder)
    # TODO: check if we should also transfer decoys
    msms_mq = msms_mq[msms_mq['Reverse'] != '+']

    logger.info(f'Reading in MaxQuant evidence.txt file and filtering out decoy hits')
    evidence_mq = mq.read_evidence_txt(mq_txt_folder)
    evidence_mq = evidence_mq[evidence_mq['Reverse'] != '+']
    rawfile_metadata = mq.get_rawfile_metadata(evidence_mq)

    logger.info(f'Reading in MaxQuant allPeptides.txt file')
    allpeptides_mq = mq.read_allpeptides_txt(mq_txt_folder)


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
        msmsscans_simsi = transfer.transfer(annotated_clusters)
        simsi_output.export_msmsscans(msmsscans_simsi, output_folder, pval)
        logger.info(f'Finished identity transfer.')

        logger.info(f'Building SIMSI-Transfer msms.txt file for {pval}.')
        msms_simsi = simsi_output.remove_unidentified_scans(msmsscans_simsi)
        simsi_output.export_msms(msms_simsi, output_folder, pval)
        logger.info(f'Finished SIMSI-Transfer msms.txt assembly.')

        statistics[pval] = simsi_output.count_clustering_parameters(msms_simsi)

        logger.info(f'Starting SIMSI-Transfer evidence.txt building for {pval}.')
        evidence_simsi = evidence.build_evidence(msms_simsi, evidence_mq, allpeptides_mq, tmt)
        simsi_output.export_simsi_evidence_file(evidence_simsi, output_folder, pval)
        logger.info(f'Finished SIMSI-Transfer evidence.txt building.')
        logger.info('')

    endtime = datetime.now()
    logger.info(f'Successfully finished transfers for all stringencies.')
    logger.info('')
    logger.info(f"SIMSI-Transfer finished in {(endtime - starttime).total_seconds()} seconds (wall clock).")


if __name__ == '__main__':
    main(sys.argv[1:])
