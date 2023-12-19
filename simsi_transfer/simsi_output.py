import os
import logging
from pathlib import Path

import pandas as pd

from simsi_transfer.merging_functions import merge_with_msmsscanstxt, merge_with_summarytxt, merge_with_msmstxt

logger = logging.getLogger(__name__)


def export_annotated_clusters(annotated_clusters, mainpath, pval):
    logger.info(f'Writing {pval}_annotated_clusters.txt.')
    export_csv(annotated_clusters, 'annotated_clusters', mainpath, pval)


def export_msmsscans(msmsscans_simsi, mainpath, pval):
    logger.info(f'Writing {pval}_msmsScans.txt for {pval}.')
    export_csv(msmsscans_simsi, 'msmsScans', mainpath, pval, sort_columns=['Raw file', 'scanID'])


def export_msms(msms_simsi, mainpath, pval):
    logger.info(f'Writing {pval}_msms.txt for {pval}.')
    export_csv(msms_simsi, 'msms', mainpath, pval, sort_columns=['Sequence', 'Modified sequence'])


def export_simsi_evidence_file(evidence_file, mainpath, pval):
    logger.info(f'Writing {pval}_evidence.txt for {pval}.')
    export_csv(evidence_file, 'evidence', mainpath, pval, sort_columns=['Sequence', 'Modified sequence'])


def export_csv(df: pd.DataFrame, filename: str, mainpath, pval, sort_columns=None):
    sumpath = mainpath / Path('summaries')
    if not os.path.exists(sumpath):
        os.makedirs(sumpath)
    pval_path = sumpath / Path(pval)
    if not os.path.exists(pval_path):
        os.makedirs(pval_path)
    if sort_columns:
        df.sort_values(by=sort_columns)
    path = pval_path / Path(f'{pval}_{filename}.txt')
    df.to_csv(path, sep='\t', index=False, na_rep='NaN')


def count_clustering_parameters(summary, rawtrans=False):
    """
    Counts MICs, tIDs, dIDs, optionally IDs with lost phospho location, and clusters. Requires flagged MICs and tIDs.
    
    MICs = multiply identified clusters, i.e. ambiguous clusters with multiple peptide sequence
    tIDs = transferred identifications (by SIMSI-Transfer)
    dIDs = direct identifications (by MaxQuant)
    
    :param summary: summary_extended input dataframe
    :param rawtrans: Flag for added lost phospho localization counting
    :return: Dictionary of counted values
    """
    ids = len(summary)
    dids = len(summary[summary['identification'] == 'd'])
    tids = len(summary[summary['identification'] == 't'])
    
    lostphos = None
    if rawtrans:
        lostphos = len(
            summary[
                (summary['identification'] == 't') & (summary['Modified sequence'] != summary['Modified sequence'])])
    totclus = max(summary['clusterID'])
    mulclus = max(summary[summary['clusterID'].duplicated(keep=False)]['clusterID'])
    pho_mics = summary[summary['mod_ambiguous'] == 1]['clusterID'].nunique(dropna=True)
    raw_mics = summary[summary['raw_ambiguous'] == 1]['clusterID'].nunique(dropna=True)
    
    logger.info(f'Identified spectra: {str(ids)}')
    logger.info(f'- MaxQuant IDs: {str(dids)}')
    logger.info(f'- SIMSI IDs: {str(tids)}')
    if rawtrans:
        logger.info(f'  - PTM-Isomeric: {str(lostphos)}')
    logger.info(f'Clusters: {str(totclus)}')
    logger.info(f'- size > 1: {str(mulclus)}')
    logger.info(f'- PTM-isomeric: {str(pho_mics)}')
    logger.info(f'- Ambiguous: {str(raw_mics)}')
    
    if rawtrans:
        return {'scans': ids, 'dids': dids, 'tids': tids, 'lostphos': lostphos, 'totclus': totclus,
                'mulclus': mulclus, 'pho_mics': pho_mics, 'raw_mics': raw_mics}
    else:
        return {'scans': ids, 'dids': dids, 'tids': tids, 'totclus': totclus, 'mulclus': mulclus,
                'pho_mics': pho_mics, 'raw_mics': raw_mics}


def remove_unidentified_scans(summary):
    summary = summary.loc[~summary['identification'].isna()]
    summary.insert(0, 'summary_ID', range(len(summary)))
    return summary


def annotate_clusters(msmsscansdf, msmsdf, rawfile_metadata, cluster_results):
    """
    Merges msmsscans.txt, msms.txt, and MaRaCluster clusters to generate summary file
    :param msmsscansdf:
    :param msmsdf:
    :param rawfile_metadata:
    :param clusterfile:
    :return:
    """
    logger.info(f'Annotating clusters with search engine results.')
    summary = merge_with_msmsscanstxt(cluster_results, msmsscansdf)
    summary = merge_with_summarytxt(summary, rawfile_metadata)
    summary = merge_with_msmstxt(summary, msmsdf)
    return summary


if __name__ == '__main__':
    raise NotImplementedError('Do not run this script.')
