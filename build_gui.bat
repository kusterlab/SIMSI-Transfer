:::: pyinstaller gui.py --noconfirm --onedir --name="SIMSI-Transfer"
pyinstaller --upx-dir ../../Downloads/upx-3.96-win64/upx-3.96-win64 --noconfirm SIMSI-Transfer.spec && python create_gui_zip.py
