import os
from pathlib import Path
from os import listdir
import shutil

simsi_dist_dir = Path.cwd() / 'dist' / 'SIMSI-Transfer'

lib_dir = simsi_dist_dir / 'lib'
lib_dir.mkdir(parents=True, exist_ok=True)

print("Moving libraries to lib directory")
for f in simsi_dist_dir.glob('*'):
    if f.name.endswith(".egg-info'") or f.name in ["pytz", "matplotlib", "sqlalchemy", "tcl8", "PIL", "greenlet", "certifi"]:
        shutil.rmtree(f)
    elif f.is_file() and f.name not in ['base_library.zip', 'python38.dll', 'python39.dll', 'SIMSI-Transfer.exe', 'pyproject.toml']:
        f.rename(lib_dir / f.name)

print("Creating zip archive")
shutil.make_archive(Path.cwd() / 'dist' / 'SIMSI-Transfer_GUI_windows', 'zip', simsi_dist_dir)
print(os.listdir())
print("Created zip file: " + str(Path.cwd() / 'dist' / 'SIMSI-Transfer_GUI_windows.zip'))