:::: pyinstaller gui.py --noconfirm --onedir --name="SIMSI-Transfer"
pyinstaller --upx-dir C:\UPX\upx-3.96\ --noconfirm ..\SIMSI-Transfer.spec && python create_gui_zip.py
