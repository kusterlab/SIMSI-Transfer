import unittest
import simsi_transfer.maracluster as cluster

class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertEqual(True, False)  # add assertion here


def test_get_file_name():
    assert cluster.get_file_name('/this/is/a/file.mzML') == 'file'


if __name__ == '__main__':
    unittest.main()
