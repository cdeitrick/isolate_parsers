import unittest
from pathlib import Path

from breseqparser.isolate_parser import get_sample_name, parse_breseq_isolate


class TestIsolateParserFunctions(unittest.TestCase):
	def setUp(self):
		self.path = Path("Clonal_Output")

	def test_get_sample_name(self):
		test = get_sample_name(Path("Clonal_Output"))
		self.assertEqual("Clonal_Output", test)

		test = get_sample_name(Path("Clonal_Output") / "breseq output")
		self.assertEqual("Clonal_Output", test)

	def test_parse_breseq_isolate(self):
		pass

if __name__ == "__main__":
	unittest.main()
