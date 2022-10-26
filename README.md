# SIMSI-Transfer

[![PyPI version](https://img.shields.io/pypi/v/simsi_transfer.svg?logo=pypi&logoColor=FFE873)](https://pypi.org/project/simsi_transfer/)
[![Supported Python versions](https://img.shields.io/pypi/pyversions/simsi_transfer.svg?logo=python&logoColor=FFE873)](https://pypi.org/project/simsi_transfer/)
[![PyPI downloads](https://img.shields.io/pypi/dm/simsi_transfer.svg)](https://pypistats.org/packages/simsi_transfer)

Transferring identifications using MS2 spectrum clustering with MaxQuant search results.

Hamood, F., Bayer, F. P., Wilhelm, M., Kuster, B., & The, M. (2022). _[SIMSI-Transfer: Software-assisted reduction of missing values in phosphoproteomic and proteomic isobaric labeling data using tandem mass spectrum clustering.](https://www.sciencedirect.com/science/article/pii/S1535947622000469)_ Molecular & Cellular Proteomics, 100238.

## Test dataset

For testing SIMSI-Transfer after installation, we recommend downloading the TMT11 MS2 raw files from this publication:
Thompson, A., Wölmer, N., Koncarevic, S., Selzer, S. et al., _[TMTpro: Design, Synthesis, and Initial Evaluation of a Proline-Based Isobaric 16-Plex Tandem Mass Tag Reagent Set.](https://pubs.acs.org/doi/abs/10.1021/acs.analchem.9b04474)_ Analytical Chemistry 2019, 91, 15941–15950. doi:10.1021/acs.analchem.9b04474

PRIDE link: https://www.ebi.ac.uk/pride/archive/projects/PXD014750

Raw files for TMT-MS2:
- 19070-001.raw
- 19070-002.raw
- 19070-003.raw
- 19070-006.raw
- 19070-007.raw
- 19070-008.raw

The MaxQuant results needed as input to SIMSI-Transfer can be downloaded from Zenodo: 
- [10.5281/zenodo.6365902](https://zenodo.org/record/6365902)

For reference, the original SIMSI-Transfer results (v0.1.0) for this dataset can also be downloaded from Zenodo:
- [10.5281/zenodo.6365638](https://zenodo.org/record/6365638)

## Running SIMSI-Transfer using the GUI

On Windows, you can download the `SIMSI-Transfer_GUI_windows.zip` from the latest release, unzip it and open `SIMSI-Transfer.exe` to start the GUI (no installation necessary).

Alternatively, on all platforms, first install SIMSI-Transfer as explained below. Then install `PyQt5` (`pip install PyQt5`) and run:

```shell
python gui.py
```

## Running SIMSI-Transfer from the command line

First install SIMSI-Transfer **as explained below**, then run SIMSI-Transfer:

```shell
python -m simsi_transfer --mq_txt_folder </path/to/txt/folder> --raw_folder </path/to/raw/folder> --output_folder </path/to/output/folder>
```

Alternative command for MS3 acquisition, using the TMT correction factor file exported from MaxQuant:

```shell
python -m simsi_transfer --mq_txt_folder </path/to/txt/folder> --raw_folder </path/to/raw/folder> --output_folder </path/to/output/folder> --tmt_ms_level ms3 --tmt_requantify --tmt_reporter_correction_file </path/to/correction/factor/file.txt>
```

Alternative command using the meta input file for MS3 acquisition, with filtered decoys:

```shell
python -m simsi_transfer --meta_input_file </path/to/meta/file> --output_folder </path/to/output/folder> --tmt_ms_level ms3 --tmt_requantify --filter_decoys
```

A list of all possible arguments is displayed using the help argument:
```shell
python -m simsi_transfer --help
```

## Installation

SIMSI-Transfer is available on PyPI and can be installed with `pip`:

```shell
pip install simsi-transfer
```

Alternatively, you can install SIMSI-Transfer after cloning from this repository:

```shell
git clone https://github.com/kusterlab/SIMSI-Transfer.git
pip install .
```
