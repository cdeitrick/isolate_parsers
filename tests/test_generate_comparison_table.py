import unittest
from math import nan

from file_generators.generate_comparison_table import *
from file_generators.generate_comparison_table import _calculate_average_value

test_snp_group = pandas.DataFrame(
	[
		{
			'alt':              'A',
			'aminoAlt':         'I',
			'aminoRef':         'M',
			'annotation':       'M350I (ATG-ATA) ',
			'codonAlt':         'ATA',
			'codonRef':         'ATG',
			'description':      'putative GGDEF domain signaling protein',
			'evidence':         'RA',
			'gene':             'PFLU3571 -',
			'index':            4,
			'locusTag':         'PFLU3571',
			'mutation':         'G-A',
			'mutationCategory': 'snp_nonsynonymous',
			'position':         3959631,
			'ref':              'G',
			'sampleId':         'S2_58BA',
			'sampleName':       'S2_58BA',
			'seq id':           'NC_012660'
		},
		{
			'alt':              'A',
			'aminoAlt':         'I',
			'aminoRef':         'M',
			'annotation':       'M350I (ATG-ATA) ',
			'codonAlt':         'ATA',
			'codonRef':         'ATG',
			'description':      'putative GGDEF domain signaling protein',
			'evidence':         'RA',
			'gene':             'PFLU3571 -',
			'index':            4,
			'locusTag':         'PFLU3571',
			'mutation':         'G-A',
			'mutationCategory': 'snp_nonsynonymous',
			'position':         3959631,
			'ref':              'G',
			'sampleId':         'S1_58BA',
			'sampleName':       'S1_58BA',
			'seq id':           'NC_012660'
		}
	]
)

test_deletion_group = pandas.DataFrame([
	{
		'alt':              'G',
		'aminoAlt':         nan,
		'aminoRef':         nan,
		'annotation':       nan,
		'codonAlt':         nan,
		'codonRef':         nan,
		'description':      'large deletion',
		'evidence':         'MC',
		'gene':             'parA-pQBR0478',
		'index':            0,
		'locusTag':         '[pQBR0001]–[pQBR0478]',
		'mutation':         'D425,094 bp',
		'mutationCategory': 'large_deletion',
		'position':         1,
		'ref':              'GATCTGCATCTGCAAAGGAGCCAGTCATGTCACCGAAGAAACGCACCCACAAGCCTGCGAAGGTCATTGTCATCGAGAACCAGAAGGGCGGCGTTGGCAA',
		'sampleId':         'SC_58SM',
		'sampleName':       'SC_58SM',
		'seq id':           'NC_009444'
	},
	{
		'alt':              'G',
		'aminoAlt':         nan,
		'aminoRef':         nan,
		'annotation':       nan,
		'codonAlt':         nan,
		'codonRef':         nan,
		'description':      'large deletion',
		'evidence':         'MC',
		'gene':             'parA-pQBR0478',
		'index':            0,
		'locusTag':         '[pQBR0001]–[pQBR0478]',
		'mutation':         'D425,094 bp',
		'mutationCategory': 'large_deletion',
		'position':         1,
		'ref':              'GATCTGCATCTGCAAAGGAGCCAGTCATGTCACCGAAGAAACGCACCCACAAGCCTGCGAAGGTCATTGTCATCGAGAACCAGAAGGGCGGCGTTGGCAA',
		'sampleId':         'S1_58BA',
		'sampleName':       'S1_58BA',
		'seq id':           'NC_009444'
	}])


class TestGenerateComparisonTable(unittest.TestCase):
	def setUp(self):
		self.maxDiff = 1000

	def test_calculate_average_value(self):
		values = ["12", "13.4", '33']

		self.assertEqual(58.4 / 3, _calculate_average_value(values))

		self.assertEqual(0.0, _calculate_average_value([]))

		self.assertEqual(0.0, _calculate_average_value(["", "", ""]))

	def test_parse_mutation_group_snp(self):
		unique_samples = ["S2_58BA", "SC_58SM", "S1_58BA"]
		group = test_snp_group
		ref_col = "ref"
		alt_col = 'alt'

		expected_output = {
			"seq id":              "NC_012660",
			"position":            3959631,
			"description":         'putative GGDEF domain signaling protein',
			"annotation":          'M350I (ATG-ATA) ',
			"locusTag":            'PFLU3571',
			"gene":                'PFLU3571 -',
			"mutationCategory":    "snp_nonsynonymous",
			"ref":                 'G',
			"S2_58BA":             "A",
			"S1_58BA":             "A",
			"SC_58SM":             "G",
			"presentInAllSamples": False,
			"presentIn":           2
		}

		output = parse_mutation_group(group, unique_samples, ref_col, alt_col)

		self.assertDictEqual(expected_output, output)

	def test_parse_mutation_group_large_deletion(self):
		unique_samples = ["S2_58BA", "SC_58SM", "S1_58BA"]
		group = test_deletion_group
		ref_col = "ref"
		alt_col = 'alt'
		ref_seq = 'GATCTGCATCTGCAAAGGAGCCAGTCATGTCACCGAAGAAA' \
				  'CGCACCCACAAGCCTGCGAAGGTCATTGTCATCGAGAACCAGAAGGGCGGCGTTGGCAA'
		expected_output = {
			"seq id":              "NC_009444",
			"position":            1,
			"description":         'large deletion',
			"annotation":          "D425,094 bp",
			"locusTag":            '[pQBR0001]–[pQBR0478]',
			"gene":                'parA-pQBR0478',
			"mutationCategory":    "large_deletion",
			"ref":                 ref_seq,
			"S2_58BA":             ref_seq,
			"S1_58BA":             "G",
			"SC_58SM":             "G",
			'presentInAllSamples': False,
			'presentIn':           2
		}

		output = parse_mutation_group(group, unique_samples, ref_col, alt_col)
		self.assertDictEqual(expected_output, output)

	def test_validate_mutation_group(self):
		pass


if __name__ == "__main__":
	pass
