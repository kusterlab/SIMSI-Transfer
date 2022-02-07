import sys
import os
import warnings
from pathlib import Path
from datetime import datetime

import pandas as pd

from .IO_functions import export_summary_file, open_msms_txt, open_msmsscans_txt, open_maracluster_clusters, \
    user_input, parse_args, open_summary_txt, open_evidence_txt, export_simsi_evidence_file
from .processing_functions import generate_summary_file, flag_ambiguous_clusters, transfer, \
    count_clustering_parameters, count_phos, build_evidence, remove_unidentified_scans

from . import thermo_raw as raw
from . import maracluster as cluster
#from . import tmt_extractor

"""
Needed:
    mainpath/
        combined/txt/.txt
Generated:
    mainpath/
        summaries/
            <pval>/<pval>summary.csv
     
            
"""


# python -m simsi_transfer.main --mq_txt_folder /media/kusterlab/internal_projects/active/Clustering_Transfers/Cluster_Tester/raw/combined/txt/ --raw_folder /media/kusterlab/internal_projects/active/Clustering_Transfers/Cluster_Tester/raw/ --output_folder /media/kusterlab/users_files/Firas_Hamood/PhD/Cluster_Tester/simsi_output
def main(argv):
    if len(argv) > 0:
        mq_txt_folder, raw_folder, pvals, output_folder, num_threads = parse_args(argv)
    else:
        mq_txt_folder, raw_folder, pvals, output_folder, num_threads = user_input()

    starttime = datetime.now()

    print(mq_txt_folder)
    print(raw_folder)
    print(pvals)
    print(output_folder)
    print(num_threads)

    mq_txt_folder, raw_folder, output_folder = Path(mq_txt_folder), Path(raw_folder), Path(output_folder)

    print(mq_txt_folder)
    print(raw_folder)
    print(pvals)
    print(output_folder)
    print(num_threads)

    print(f'{datetime.now()}:\tStarting SIMSI-Transfer')
    print()

    if not output_folder.is_dir():
        output_folder.mkdir(parents=True)

    print(f'{datetime.now()}:\tConverting .raw files')
    mzml_folder = output_folder / Path('mzML')
    mzml_files = raw.convert_raw_mzml_batch(raw_folder, mzml_folder, num_threads)

    #print(f'{datetime.now()}:\tExtract reporter ion intensities from .mzML files')
    #extracted_folder = output_folder / Path('extracted')
    #tmt_extractor.extract_tmt_reporters(mzml_files, extracted_folder, num_threads)

    print(f'{datetime.now()}:\tClustering .mzML files')
    maracluster_folder = output_folder / Path('maracluster_output')
    cluster.cluster_mzml_files(mzml_files, pvals, maracluster_folder, num_threads)

    ####################################################################################################################

    print(f'{datetime.now()}:\tOpening MaxQuant msmsscans.txt file and clearing columns')
    msmsscanstxt, tmt = open_msmsscans_txt(mq_txt_folder)

    print(f'{datetime.now()}:\tOpening MaxQuant msms.txt file and filtering out decoy hits')
    msmstxt = open_msms_txt(mq_txt_folder)
    # TODO: clean
    #msmstxt = msmstxt[msmstxt['Reverse'] != '+']

    print(f'{datetime.now()}:\tOpening MaxQuant evidence.txt file and filtering out decoy hits')
    evidencetxt = open_evidence_txt(mq_txt_folder)
    # TODO: clean
    #evidencetxt = evidencetxt[evidencetxt['Reverse'] != '+']

    print(f'{datetime.now()}:\tOpening MaxQuant summary.txt file')
    summarytxt = open_summary_txt(mq_txt_folder)
    print(summarytxt)

    statistics = dict()

    for pval in ['p' + str(-1*i) for i in pvals]:
        print()

        print(f'{datetime.now()}:\tStarting MaxQuant and MaRaCluster file merge for {pval}.')

        clusterfile = open_maracluster_clusters(maracluster_folder, pval)
        summary = generate_summary_file(msmsscanstxt, msmstxt, summarytxt, clusterfile)
        del clusterfile

        export_summary_file(summary, output_folder, pval, state='merged')
        print(f'{datetime.now()}:\tFinished file merge.')

        print(f'{datetime.now()}:\tStarting cluster-based identity transfer for {pval}.')
        summary['phosphogroups'] = summary['Modified sequence'].apply(count_phos)
        summary = flag_ambiguous_clusters(summary)
        summary = transfer(summary)
        export_summary_file(summary, output_folder, pval, state='transferred')
        print(f'{datetime.now()}:\tFinished identity transfer.')

        print(f'{datetime.now()}:\tBuilding msms.csv file for {pval}.')
        summary = remove_unidentified_scans(summary)
        export_summary_file(summary, output_folder, pval, state='filtered')
        print(f'{datetime.now()}:\tFinished msms.csv assembly.')

        statistics[pval] = count_clustering_parameters(summary)

        print(f'{datetime.now()}:\tStarting evidence_transferred.txt building for {pval}.')
        evidence = build_evidence(summary, evidencetxt, tmt)
        export_simsi_evidence_file(evidence, output_folder, pval)
        print(f'{datetime.now()}:\tFinished evidence_transferred.txt building.')
        print()

    endtime = datetime.now()
    print(f'{datetime.now()}:\tSuccessfully finished file processing for all stringencies.')
    print()
    print(f'{endtime - starttime} needed for full run.')
    print(f'Terminating SIMSI-Transfer...')


if __name__ == '__main__':
    main(sys.argv[1:])
