import os
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def export_summary_file(summary_file, mainpath, pval, state):
    """
    Generate and return summary dataframes, a merge of msms.txt identifications and MaRaCluster clustering results with
    one line for each ms2 scan generated during acquisition
    :param summary_file: summary of msms.txt, msmsscans.txt, and MaRaCluster cluster.csv files
    :param mainpath: Path to main processing folder
    :param pval: MaRaCluster stringency threshold; a higher value (e.g. p30) is more stringent than a lower one (e.g.
    p15), but results in less transfers
    :param state: set to 'merged', 'transferred', or 'filtered' for merged files, msmsScans, and msms files respectively
    """
    sumpath = mainpath / Path('summaries')
    if not os.path.exists(sumpath):
        os.makedirs(sumpath)
    pval_path = sumpath / Path(pval)
    if not os.path.exists(pval_path):
        os.makedirs(pval_path)
    if state == 'merged':
        summary_file.to_csv(pval_path / Path(f'{pval}_summary.csv'), sep='\t', index=False, na_rep='NaN')
    elif state == 'transferred':
        summary_file = summary_file.sort_values(by=['Raw file', 'scanID'])
        summary_file.to_csv(pval_path / Path(f'{pval}_msmsScans.csv'), sep='\t', index=False, na_rep='NaN')
    elif state == 'filtered':
        summary_file = summary_file.sort_values(by=['Sequence', 'Modified sequence'])
        summary_file.to_csv(pval_path / Path(f'{pval}_msms.csv'), sep='\t', index=False, na_rep='NaN')
    else:
        raise ValueError("Incorrect state given, please set state parameter to 'merged', 'transferred', or 'filtered'")


def export_simsi_evidence_file(evidence_file, mainpath, pval):
    sumpath = mainpath / Path('summaries')
    pval_path = sumpath / Path(pval)
    evidence_file.to_csv(pval_path / Path(f'{pval}_evidence.csv'), sep='\t', index=False, na_rep='NaN')





if __name__ == '__main__':
    raise NotImplementedError('Do not run this script.')


def count_clustering_parameters(summary, rawtrans=False):
    """
    Counts MICs, tIDs, dIDs, optionally IDs with lost phospho location, and clusters. Requires flagged MICs and tIDs.
    :param summary: summary_extended input dataframe
    :param rawtrans: Flag for added lost phospho localization counting
    :return: Dictionary of counted values
    """
    scans = len(summary)
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
    logger.info(f'Scans: {str(scans)}')
    logger.info(f'MaxQuant IDs: {str(dids)}')
    logger.info(f'SIMSI IDs: {str(tids)}')
    if rawtrans:
        logger.info(f'Isomeric SIMSI IDs: {str(lostphos)}')
    logger.info(f'All clusters: {str(totclus)}')
    logger.info(f'Clusters > 1: {str(mulclus)}')
    logger.info(f'PTM-isomeric clusters : {str(pho_mics)}')
    logger.info(f'Ambiguous clusters: {str(raw_mics)}')
    if rawtrans:
        return {'scans': scans, 'dids': dids, 'tids': tids, 'lostphos': lostphos, 'totclus': totclus,
                'mulclus': mulclus, 'pho_mics': pho_mics, 'raw_mics': raw_mics}
    else:
        return {'scans': scans, 'dids': dids, 'tids': tids, 'totclus': totclus, 'mulclus': mulclus,
                'pho_mics': pho_mics, 'raw_mics': raw_mics}


def remove_unidentified_scans(summary):
    summary = summary.loc[~summary['identification'].isna()]
    summary.insert(0, 'summary_ID', range(len(summary)))
    return summary