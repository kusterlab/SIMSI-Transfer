# clustering_transfer_tool

Tool for increasing PSM gain from MaxQuant output file. Requires "Experiment" set in MaxQuant!


Stuff left to do:
- Argparse instead of input()
- Implement masking analysis for FDR estimation (?)
- Add mzML conversion step to pipeline (?)
- SIMSI-GUI

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

3. We need to install scipy through `pip`, because the conda version depends on the MKL library again. Here, we also need an older version of scipy (see [here](https://stackoverflow.com/questions/62581504/why-do-i-have-modulenotfounderror-no-module-named-scipy-special-cython-specia))

```
pip install scipy==1.4.1
```

## Building a self-contained executable

Use the `build.bat` script to create a self-contained executable. Note that we need some adjustments to the command line for [hidden imports](https://stackoverflow.com/questions/57108026/pyinstaller-modulenotfounderror-no-module-named-sklearn-utils-cython-blas) for `sklearn`.


## Reducing size of the executable

Using UPX (https://upx.github.io/) reduces the DLL file sizes:

Compress everything except these DLLs:

* msvcp140.dll
* VCRUNTIME140.dll
* PyQt5\Qt\plugins\platforms\windows.dll
* scipy\_lib\_uarray\_uarray.cp36-win_amd64.pyd

```
call upx.bat
```

We can exclude some other modules/DLLs to reduce the file size even further:

* OpenGL and EGL can be excluded if we remove the DLL loading in `C:\ProgramData\Anaconda3\envs\simsi_transfer_gui\Lib\site-packages\mPyInstaller\utils\hooks\qt.py` (lines 609-622)
