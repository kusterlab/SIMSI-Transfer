# clustering_transfer_tool

Tool for increasing PSM gain from MaxQuant output file. Requires "Experiment" set in MaxQuant!

Example usage:
```shell
python -m simsi_transfer.main --mq_txt_folder </path/to/txt/folder> --raw_folder </path/to/raw/folder> --output_folder </path/to/output/folder>
```


Stuff left to do:
- Implement masking analysis for FDR estimation (?)

# Building the GUI on Windows

## Create a conda environment

Try importing `conda_environment.yml` in the Anaconda environment tab.

If that does not work, try the following:

1. Set up a new environment, either through the Anaconda UI, or by running the following on the command line:

```
conda create -n simsi_transfer_gui python=3.6
activate simsi_transfer_gui
```

2. There are some caveats with installing the dependencies. We want to avoid dependence on the MKL (Math Kernel Library) package by numpy/scipy, as this blows up the size of the .exe file over 200MB (see [here](https://github.com/pyinstaller/pyinstaller/issues/2270)).

```
conda install -c conda-forge nomkl numpy pandas pyqt pyinstaller
```

## Building a self-contained executable

Use the `build_gui.bat` script to create a self-contained executable.


### Reducing size of the executable

Download UPX (https://upx.github.io/) to reduce the DLL file sizes, change the path in `build_gui.bat` to point to the UPX **folder**.
