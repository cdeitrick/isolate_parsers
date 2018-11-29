import unittest
from pathlib import Path

from breseqparser.isolate_parser import get_sample_name, parse_breseq_isolate, IsolateTableColumns

data_folder = Path(__file__).parent / 'data' / 'Clonal_Output' / 'breseq_output'

class TestIsolateParserFunctions(unittest.TestCase):
	def setUp(self):
		self.path = Path("Clonal_Output")

	def test_get_sample_name(self):
		test = get_sample_name(Path("Clonal_Output"))
		self.assertEqual("Clonal_Output", test)

		test = get_sample_name(Path("Clonal_Output") / "breseq output")
		self.assertEqual("Clonal_Output", test)

	def test_parse_breseq_isolate(self):
		breseq_folder = Path(__file__).parent / "Clonal_Output"
		variant_table, coverage_table, junction_table = parse_breseq_isolate(breseq_folder, isolate_id = 'testIsolate')

		self.assertListEqual(IsolateTableColumns, list(variant_table.columns))

if __name__ == "__main__":
	unittest.main()
