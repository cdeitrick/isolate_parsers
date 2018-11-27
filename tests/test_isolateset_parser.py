import unittest
from pathlib import Path
from isolateset_parser import _parse_commandline_list, _parse_sample_map
data_folder = Path(__file__).with_name('data')

class TestCommandlineParser(unittest.TestCase):
	def test_parse_sample_map(self):
		sample_map_file = data_folder / "test_sample_map.txt"

		truth_data = {
			'sample1': 'sampleA',
			'sample2': 'sampleB',
			'sample3': 'sampleC',
			'sample4': '123'
		}
		test_data = _parse_sample_map(sample_map_file)

		self.assertDictEqual(truth_data, test_data)

	def test_parse_commandline_list(self):
		teststring = "a,b,c"

		self.assertListEqual(teststring.split(','), _parse_commandline_list(teststring))


if __name__ == "__main__":
	unittest.main()