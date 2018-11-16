import unittest
from pathlib import Path

import pandas

from breseqparser.file_parsers import parse_index


class TestSetup(unittest.TestCase):
	breseq_run_folder = Path(__file__).with_name("breseq_run")
	breseq_sample_folder = breseq_run_folder / "AU0074"
	filename = breseq_sample_folder / "breseq_output" / "output" / "index.html"

	snp, coverage, junction = parse_index.parse_index_file('AU0074', filename)


class TestGetFilename(unittest.TestCase):
	def setUp(self):
		self.breseq_run_folder = Path(__file__).with_name("breseq_run")
		self.breseq_sample_folder = self.breseq_run_folder / "AU0074"
		self.filename = self.breseq_sample_folder / "breseq_output" / "output" / "index.html"

	def test_breseq_run_folder(self):
		self.assertRaises(FileNotFoundError, parse_index.get_index_filename, self.breseq_run_folder)

	def test_sample_run_folder(self):
		testpath = parse_index.get_index_filename(self.breseq_sample_folder)
		self.assertEqual(self.filename, testpath)

	def test_index_filename(self):
		testpath = parse_index.get_index_filename(self.filename)
		self.assertEqual(self.filename, testpath)


class TestIndexTable(TestSetup):
	def test_snp_is_dataframe(self):
		self.assertIsInstance(self.snp, pandas.DataFrame)
		self.assertEqual(1070, len(self.snp))

	def test_coverage_is_dataframe(self):
		self.assertIsInstance(self.coverage, pandas.DataFrame)
		self.assertEqual(18, len(self.coverage))

	def test_junction_is_dataframe(self):
		self.assertIsInstance(self.junction, pandas.DataFrame)
		self.assertEqual(46, len(self.junction))

	def test_snp_has_correct_columns(self):
		expected_columns = {
			'Sample', 'evidence', 'seq id', 'position', 'mutation',
			'annotation', 'gene', 'description'
		}
		self.assertSetEqual(expected_columns, set(self.snp.columns))

	def test_coverage_has_correct_columns(self):
		expected_columns = {
			'Sample', 'seq id', 'start', 'end', 'size', '←reads', 'reads→',
			'gene', 'description', '\xa0'
		}

		self.assertSetEqual(expected_columns, set(self.coverage.columns))

	def test_junction_has_correct_columns(self):
		expected_columns = {
			'Sample', 'seq id', 'position', 'score', 'skew',
			'freq', 'annotation', 'gene', 'product', 'reads (cov)',
			'reads (cov) (single)', '0', '1'
		}
		self.assertSetEqual(expected_columns, set(self.junction.columns))


class TestSNPColumnValues(TestSetup):
	def test_sample_column(self):
		self.assertIn('Sample', self.snp.columns)
		series = self.snp['Sample']
		self.assertEqual(object, series.dtype)
		self.assertIsInstance(series.iloc[0], str)
		self.assertEqual(1, len(series.unique()))

	def test_evidence_column(self):
		self.assertIn('evidence', self.snp.columns)
		series = self.snp['evidence']

		self.assertIsInstance(series.iloc[0], str)
		unique = set(series.unique().tolist())
		self.assertSetEqual({'RA', 'JC', 'MC JC'}, unique)

	def test_seq_column(self):
		self.assertIn('seq id', self.snp.columns)
		series = self.snp['seq id']
		self.assertEqual(object, series.dtype)
		self.assertIsInstance(series.iloc[0], str)
		self.assertGreaterEqual(len(series.unique()), 1)

	def test_position_column(self):
		self.assertIn('position', self.snp.columns)
		series = self.snp['position']
		self.assertEqual(int, series.dtype)

	def test_mutation_column(self):
		self.assertIn('mutation', self.snp.columns)
		series = self.snp['mutation']
		self.assertIsInstance(series.iloc[0], str)


if __name__ == "__main__":
	unittest.main()
