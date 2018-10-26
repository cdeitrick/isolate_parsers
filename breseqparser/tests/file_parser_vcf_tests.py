import unittest
from pathlib import Path
from breseqparser.file_parsers import index_file_parser, vcf_file_parser
class TestGetFilename(unittest.TestCase):
	def setUp(self):
		self.truth = Path(__file__).parent / "AU0074" / "breseq_output" / "output" / "index.html"
	def test_single_run(self):
		folder = Path(__file__).parent
		testpath = index_file_parser.get_index_filename(folder)
		self.assertEqual(self.truth, testpath)

	def test_direct_link(self):
		testpath = index_file_parser.get_index_filename(self.truth)
		self.assertEqual(self.truth, testpath)

	def test_link_to_run(self):
		link = Path(__file__).parent / "AU0074"
		testpath = index_file_parser.get_index_filename(link)
		self.assertEqual(self.truth, testpath)

if __name__ == "__main__":
	unittest.main()

