import unittest
from pathlib import Path

import pandas
import itertools
from breseqparser.file_parsers import parse_vcf


class TestSetup(unittest.TestCase):
	def setUp(self):
		self.breseq_run_folder = Path(__file__).with_name("breseq_run")
		self.breseq_sample_folder = self.breseq_run_folder / "AU0074"
		self.filename = self.breseq_sample_folder / "breseq_output" / "data" / "output.vcf"
		self.output = parse_vcf.parse_vcf(self.filename)
		self.output.reset_index(inplace = True)


class TestGetVCFFilename(TestSetup):
	def test_breseq_run_folder(self):
		self.assertRaises(FileNotFoundError, parse_vcf.get_vcf_filename, self.breseq_run_folder)

	def test_sample_run_folder(self):
		testpath = parse_vcf.get_vcf_filename(self.breseq_sample_folder)
		self.assertEqual(self.filename, testpath)

	def test_index_filename(self):
		testpath = parse_vcf.get_vcf_filename(self.filename)
		self.assertEqual(self.filename, testpath)


class TestVCFTable(TestSetup):
	def test_is_dataframe(self):
		self.assertIsInstance(self.output, pandas.DataFrame)
		self.assertEqual(len(self.output), 1070)

	def test_has_correct_columns(self):
		expected_columns = {
			'sampleName', 'seq id', 'position', 'alt', 'ref',
			'quality', 'readDepth', 'variantType'
		}
		self.assertSetEqual(expected_columns, set(self.output.columns))


class TestVCFTableColumnValues(TestSetup):
	def test_samplename(self):
		self.assertIn('sampleName', self.output.columns)
		series = self.output['sampleName']
		self.assertEqual(object, series.dtype)
		self.assertIsInstance(series.iloc[0], str)
		self.assertEqual(1, len(series.unique()))

	def test_seq_id(self):
		self.assertIn('seq id', self.output.columns)
		series = self.output['seq id']
		self.assertEqual(object, series.dtype)
		self.assertIsInstance(series.iloc[0], str)
		self.assertGreaterEqual(len(series.unique()), 1)

	def test_ref_column(self):
		self.assertIn('ref', self.output.columns)
		series = self.output['ref']
		self.assertIsInstance(series.iloc[0], str)
		unique = set(itertools.chain.from_iterable(series.unique().tolist()))
		self.assertSetEqual(set(), unique-{'A', 'C', 'G', 'T', 'N', 'U', '.'})

	def test_alt_column(self):
		self.assertIn('alt', self.output.columns)
		series = self.output['alt']
		self.assertIsInstance(series.iloc[0], str)
		unique = set(itertools.chain.from_iterable(series.unique().tolist()))
		self.assertSetEqual(set(), unique-{'A', 'C', 'G', 'T', 'N', 'U', '.'})

	def test_quality_column(self):
		self.assertIn('quality', self.output.columns)
		series = self.output['quality']

		self.assertEqual(float, series.dtype)
		self.assertIsInstance(series.iloc[0], float)

	def test_readdepth_column(self):
		self.assertIn('readDepth', self.output.columns)
		series = self.output['readDepth']

		self.assertEqual(float, series.dtype)
		self.assertIsInstance(series.iloc[0], float)


if __name__ == "__main__":
	unittest.main()
