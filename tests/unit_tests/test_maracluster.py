import simsi_transfer.maracluster as cluster


def test_get_file_name():
    assert cluster.get_file_name('/this/is/a/file.mzML') == 'file'

