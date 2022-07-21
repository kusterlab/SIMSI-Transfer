import simsi_transfer.thermo_raw as raw


def test_apply_and_flatten():
    raw_folders = ['folder_1', 'folder_2', 'folder_3']
    def list_files(folder):
        return [f'{folder}/file_1.txt', f'{folder}/file_1.txt', f'{folder}/file_2.txt']
    
    raw_files = raw.apply_and_flatten(raw_folders, list_files)
    assert raw_files == ['folder_1/file_1.txt', 'folder_1/file_2.txt', 'folder_2/file_1.txt', 'folder_2/file_2.txt', 'folder_3/file_1.txt', 'folder_3/file_2.txt']
        
