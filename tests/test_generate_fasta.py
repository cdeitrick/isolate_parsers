import unittest

from file_generators.generate_fasta import *
# Need to also import private functions
from file_generators.generate_fasta import _get_relevant_columns, _parse_sample_group, _validate_variant_table

data_folder = Path(__file__).parent / 'data' / 'Clonal_Output' / 'breseq_output'
gd_filename = data_folder / 'output' / 'evidence' / 'annotated.gd'
test_table = [
	{'ref': 'G', 'alt': 'A', 'sampleName': 'S2_58BA', 'seq id': 'NC_012660', 'position': 3959631},
	{'ref': 'A', 'alt': 'C', 'sampleName': 'SC_58SM', 'seq id': 'NC_012660', 'position': 4173231},
	{'ref': 'G', 'alt': 'A', 'sampleName': 'S1_58BA', 'seq id': 'NC_012660', 'position': 3959631}
]


class TestFastaGenerator(unittest.TestCase):
	def test_get_relevant_columns(self):
		ref_base, alt_base = _get_relevant_columns('base')
		self.assertEqual('ref', ref_base)
		self.assertEqual('alt', alt_base)

		ref_codon, alt_codon = _get_relevant_columns('codon')
		self.assertEqual('codonRef', ref_codon)
		self.assertEqual('codonAlt', alt_codon)

		ref_amino, alt_amino = _get_relevant_columns('amino')
		self.assertEqual('aminoRef', ref_amino)
		self.assertEqual('aminoAlt', alt_amino)

	def test_generate_reference_sequence(self):
		reference_column = 'ref'
		test_df = pandas.DataFrame(test_table)
		truth_reference: pandas.Series = pandas.DataFrame(
			[
				{'seq id': 'NC_012660', 'position': 3959631, 'reference': 'G'},
				{'seq id': 'NC_012660', 'position': 4173231, 'reference': 'A'}
			]
		).set_index(keys = ['seq id', 'position'])['reference']

		test_reference = generate_reference_sequence(test_df, 'ref')

		pandas.testing.assert_series_equal(truth_reference, test_reference)

	def test_validate_variant_table(self):
		# Make sure the test table passes
		test_df = pandas.DataFrame(
			test_table + [{'ref': 'GG', 'alt': 'A', 'sampleName': 'S2_58BA', 'seq id': 'NC_012660', 'position': 3959631}]
		)
		_validate_variant_table(test_df[:-1], 'base', 'ref', 'alt')
		with self.assertRaises(ValueError):
			_validate_variant_table(test_df, 'base', 'ref', 'alt')

	def test_convert_combined_table_to_aligned_table(self):
		pass

	def test_parse_sample_group(self):
		sample_group = (
			'SC_58SM',
			pandas.DataFrame(
				[
					{'ref': 'A', 'alt': 'C', 'sampleName': 'SC_58SM', 'seq id': 'NC_012660', 'position': 4173231}
				]
			)
		)

		truth_reference: pandas.Series = pandas.DataFrame(
			[
				{'seq id': 'NC_012660', 'position': 3959631, 'reference': 'G'},
				{'seq id': 'NC_012660', 'position': 4173231, 'reference': 'A'}
			]
		).set_index(keys = ['seq id', 'position'])['reference']

		_truth_output_index = pandas.MultiIndex(
			levels = [['NC_012660'], [3959631, 4173231]],
			labels = [[0, 0], [0, 1]],
			names = ['seq id', 'position']
		)
		truth_output = pandas.Series(
			data = ["G", "C"],
			index = _truth_output_index,
			name = "SC_58SM"
		)

		output = _parse_sample_group(sample_group[0], sample_group[1], truth_reference, 'alt')

		pandas.testing.assert_series_equal(truth_output, output)


if __name__ == "__main__":
	unittest.main()
