import unittest

from breseqparser.file_parsers.parse_index import *
# Need to also import private functions
from breseqparser.file_parsers.parse_index import _extract_index_file_tables, _load_index_file

data_folder = Path(__file__).parent / 'data' / 'Clonal_Output' / 'breseq_output'
index_filename = data_folder / 'output' / 'index.html'
snp_table_filename = data_folder / "snp_table.pkl"
coverage_table_filename = data_folder / "coverage_table.pkl"
junction_table_filename = data_folder / "junction_table.pkl"


class TestIndexParser(unittest.TestCase):
	def setUp(self):
		# Preload the Beautifulsoup object to reduce runtimes.
		self.soup = _load_index_file(index_filename)
		self.snp_table = pandas.read_pickle(snp_table_filename)
		self.coverage_table = pandas.read_pickle(coverage_table_filename)
		self.junction_table = pandas.read_pickle(junction_table_filename)

	def test_convert_to_number(self):
		self.assertEqual(123456, to_number('123,456'))
		self.assertEqual(123456, to_number('123456'))

	def test_load_index_file(self):
		soup = _load_index_file(index_filename)
		self.assertIsInstance(soup, BeautifulSoup)

	def test_extract_index_file_tables(self):
		# Omit `seq id` since there is only one sequence.
		columns = ['evidence', 'position', 'mutation', 'annotation', 'gene', 'description']

		test_headers, test_coverage, test_junction = _extract_index_file_tables(self.soup)

		self.assertListEqual(columns, test_headers)
		self.assertIsInstance(test_coverage, BeautifulSoup)
		self.assertIsInstance(test_coverage, BeautifulSoup)

	def test_index_parser_output(self):
		snp_table, coverage_table, junction_table = parse_index_file('testIsolate', index_filename)

		pandas.testing.assert_frame_equal(self.snp_table, snp_table)
		pandas.testing.assert_frame_equal(self.coverage_table, coverage_table)
		pandas.testing.assert_frame_equal(self.junction_table, junction_table)


if __name__ == "__main__":
	pass
