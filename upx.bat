cd dist\PickedGroupFDR
for /r %%e in (*.exe,*.dll,*.pyd) do ..\..\..\Downloads\upx-3.96-win64\upx-3.96-win64\upx.exe "%%e" --best 