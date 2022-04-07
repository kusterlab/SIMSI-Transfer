import os
import logging
from pathlib import Path

import pandas as pd

from .processing_functions import purge_mrc_files


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