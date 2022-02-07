:::pyinstaller main.py --noconfirm --hidden-import="sklearn.utils._cython_blas" --exclude-module matplotlib --upx-dir="C:\Users\mthe\Downloads\upx-3.96-win64\upx-3.96-win64" --upx-exclude=ucrtbase.dll --upx-exclude=vcruntime140.dll --name="SIMSI-Transfer"
pyinstaller main.py --noconfirm --hidden-import="sklearn.utils._cython_blas" --exclude-module matplotlib --name="SIMSI-Transfer"
