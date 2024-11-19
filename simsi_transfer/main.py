import sys
import os
from pathlib import Path
from datetime import datetime
import logging
import time

from . import __version__, __copyright__
from . import command_line_interface as cli
from . import maxquant as mq
from . import simsi_output
from . import thermo_raw as raw
from . import maracluster as cluster
from . import tmt_processing
from . import transfer
from . import evidence
from .utils import utils

logger = logging.getLogger(__name__)


def main(argv):
    args = cli.parse_args(argv)

    pvals = cli.parse_stringencies(args.stringencies)
    meta_input_df = cli.get_input_folders(args)
    tmt_ms_level = cli.parse_tmt_ms_level(args.tmt_ms_level)

    raw_folders = utils.convert_to_path_list(meta_input_df['raw_folder'])
    mq_txt_folders = utils.convert_to_path_list(meta_input_df['mq_txt_folder'])
    tmt_correction_files = utils.convert_to_path_list(meta_input_df['tmt_correction_file'])

    args.output_folder.mkdir(parents=True, exist_ok=True)
    args.cache_folder.mkdir(parents=True, exist_ok=True)

    module_name = ".".join(__name__.split(".")[:-1])
    file_logger = logging.FileHandler(args.output_folder / Path('SIMSI.log'))
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
    logger.info(f"Output folder = {args.output_folder}")
    logger.info(f"Cache folder = {args.cache_folder}")
    logger.info(f"Number of threads = {args.num_threads}")
    logger.info(f"TMT correction file = {tmt_correction_files}")
    logger.info(f"TMT MS level = {tmt_ms_level}")
    logger.info('')

    logger.info(f'Starting SIMSI-Transfer')
    logger.info('')

    logger.info(f'Retrieving .raw files')
    meta_input_df['raw_files'] = meta_input_df['raw_folder'].apply(raw.get_raw_files)
    raw_files, correction_factor_paths = utils.get_raw_files_and_correction_factor_paths(meta_input_df)

    raw_filenames_input = {i.stem for i in raw_files}

    logger.info(f'Converting .raw files')
    mzml_folder = args.cache_folder / Path('mzML')
    mzml_files = raw.convert_raw_mzml_batch(raw_files, mzml_folder, args.num_threads)

    cluster_result_folder = args.output_folder / Path('maracluster_output')
    if not cluster.has_previous_run(cluster_result_folder, mzml_files, pvals):
        logger.info(f'Clustering .mzML files')
        dat_files_folder = args.cache_folder / Path('dat_files')
        cluster.cluster_mzml_files(mzml_files, pvals, cluster_result_folder, dat_files_folder, args.num_threads)
    else:
        logger.info("Found previous MaRaCluster run, skipping clustering")

    plex = mq.get_plex(mq_txt_folders)

    logger.info(f'Reading in MaxQuant msmsscans.txt file')
    msmsscans_mq = utils.process_and_concat(mq_txt_folders, mq.read_msmsscans_txt, tmt_requantify=args.tmt_requantify,
                                         plex=plex)

    raw_filenames_mq = set(msmsscans_mq['Raw file'].unique())
    if raw_filenames_mq != raw_filenames_input:
        raise ValueError(
            f'The raw files listed as input and the raw files in the MaxQuant search results are not the same!')

    if args.tmt_requantify:
        logger.info(f'Extracting correct reporter ion intensities from .mzML files')
        extracted_folder = args.cache_folder / Path('extracted')
        tmt_processing.extract_tmt_reporters(mzml_files=mzml_files, output_path=extracted_folder,
                                             correction_factor_paths=correction_factor_paths, plex=plex,
                                             extraction_level=tmt_ms_level, num_threads=args.num_threads)

        corrected_tmt = tmt_processing.assemble_corrected_tmt_table(mzml_files, extracted_folder, plex)
        msmsscans_mq = tmt_processing.merge_with_corrected_tmt(msmsscans_mq, corrected_tmt)

    logger.info(f'Reading in MaxQuant msms.txt file')
    msms_mq = utils.process_and_concat(mq_txt_folders, mq.read_msms_txt)
    if args.filter_decoys:
        logger.info(f'Filtering out decoy hits')
        msms_mq = msms_mq[msms_mq['Reverse'] != '+']

    logger.info(f'Reading in MaxQuant evidence.txt file')
    evidence_mq = utils.process_and_concat(mq_txt_folders, mq.read_evidence_txt)
    if args.filter_decoys:
        logger.info(f'Filtering out decoy hits')
        evidence_mq = evidence_mq[evidence_mq['Reverse'] != '+']
    rawfile_metadata = mq.get_rawfile_metadata(evidence_mq)

    logger.info(f'Reading in MaxQuant allPeptides.txt file')
    allpeptides_mq = utils.process_and_concat(mq_txt_folders, mq.read_allpeptides_txt)
    allpeptides_mq = mq.fill_missing_min_max_scans(allpeptides_mq, msmsscans_mq)

    statistics = dict()

    for pval in ['p' + str(i) for i in pvals]:
        logger.info('')

        logger.info(f'Starting MaxQuant and MaRaCluster file merge for {pval}.')

        cluster_results = cluster.read_cluster_results(cluster_result_folder, pval)
        annotated_clusters = simsi_output.annotate_clusters(msmsscans_mq, msms_mq, rawfile_metadata, cluster_results)
        del cluster_results

        if not args.skip_annotated_clusters:
            simsi_output.export_annotated_clusters(annotated_clusters, args.output_folder, pval)
        logger.info(f'Finished file merge.')

        if args.skip_msmsscans and args.skip_msms and args.skip_evidence:
            del annotated_clusters
            continue

        logger.info(f'Starting cluster-based identity transfer for {pval}.')
        annotated_clusters = transfer.flag_ambiguous_clusters(annotated_clusters)
        msmsscans_simsi = transfer.transfer(annotated_clusters, ambiguity_decision=args.ambiguity_decision, max_pep=args.maximum_pep)
        del annotated_clusters

        if not args.skip_msmsscans:
            simsi_output.export_msmsscans(msmsscans_simsi, args.output_folder, pval)
        logger.info(f'Finished identity transfer.')

        if args.skip_msms and args.skip_evidence:
            del msmsscans_simsi
            continue

        logger.info(f'Building SIMSI-Transfer msms.txt file for {pval}.')
        msms_simsi = simsi_output.remove_unidentified_scans(msmsscans_simsi)
        del msmsscans_simsi
        
        if args.add_plotting_columns:
            raise NotImplementedError()
        
        if not args.skip_msms:
            simsi_output.export_msms(msms_simsi, args.output_folder, pval)
        logger.info(f'Finished SIMSI-Transfer msms.txt assembly.')

        statistics[pval] = simsi_output.count_clustering_parameters(msms_simsi)

        if not args.skip_evidence:
            logger.info(f'Starting SIMSI-Transfer evidence.txt building for {pval}.')
            evidence_simsi = evidence.build_evidence_grouped(msms_simsi, evidence_mq, allpeptides_mq, plex, num_threads=args.num_threads)
            simsi_output.export_simsi_evidence_file(evidence_simsi, args.output_folder, pval)
            logger.info(f'Finished SIMSI-Transfer evidence.txt building.')
            logger.info('')

        del msms_simsi, evidence_simsi

    endtime = datetime.now()
    logger.info(f'Successfully finished transfers for all stringencies.')
    logger.info('')
    logger.info(f"SIMSI-Transfer finished in {(endtime - starttime).total_seconds()} seconds (wall clock).")


if __name__ == '__main__':
    main(sys.argv[1:])
