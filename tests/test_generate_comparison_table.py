from math import nan

import pandas
import pytest

from dataio import import_table
from isolateparser.file_generators import generate_comparison_table


@pytest.fixture
def snp_group() -> pandas.DataFrame:
	test_snp_group = pandas.DataFrame(
		{
			'alt':              ['A', 'A'],
			'aminoAlt':         ['I', 'I'],
			'aminoRef':         ['M', 'M'],
			'annotation':       ['M350I (ATG-ATA) ', 'M350I (ATG-ATA) '],
			'codonAlt':         ['ATA', 'ATA'],
			'codonRef':         ['ATG', 'ATG'],
			'description':      ['putative GGDEF domain signaling protein', 'putative GGDEF domain signaling protein'],
			'evidence':         ['RA', 'RA'],
			'gene':             ['PFLU3571 -', 'PFLU3571 -'],
			'index':            [4, 4],
			'locusTag':         ['PFLU3571', 'PFLU3571'],
			'mutation':         ['G-A', 'G-A'],
			'mutationCategory': ['snp_nonsynonymous', 'snp_nonsynonymous'],
			'position':         [3959631, 3959631],
			'ref':              ['G', 'G'],
			'sampleId':         ['S2_58BA', 'S1_58BA'],
			'sampleName':       ['S2_58BA', 'S1_58BA'],
			'seq id':           ['NC_012660', 'NC_012660']
		}
	)
	return test_snp_group

@pytest.fixture
def mutation_group() -> pandas.DataFrame:
	string = """
		alt	aminoAlt	aminoRef	annotation	codonAlt	codonRef	description	evidence	gene	locusTag	mutation	mutationCategory	position	ref	sampleId	sampleName	seq id
		A	M	L	L701M (CTT-ATG)	ATG	CTT	amino acid adenylation protein	RA	BPFJKKCJ_01788 -		C-A	snp_nonsynonymous	641399	C	AU1581	E-01	NODE_2
		A	M	L	L701M (CTT-ATG)	ATG	CTT	amino acid adenylation protein	RA	BPFJKKCJ_01788 -		C-A	snp_nonsynonymous	641399	C	AU1836	E-02	NODE_2
		A	M	L	L701M (CTT-ATG)	ATG	CTT	amino acid adenylation protein	RA	BPFJKKCJ_01788 -		C-A	snp_nonsynonymous	641399	C	AU4381	E-05	NODE_2
		A	M	L	L701M (CTT-ATG)	ATG	CTT	amino acid adenylation protein	RA	BPFJKKCJ_01788 -		C-A	snp_nonsynonymous	641399	C	AU4993	E-06	NODE_2
		A	M	L	L701M (CTT-ATG)	ATG	CTT	amino acid adenylation protein	RA	BPFJKKCJ_01788 -		C-A	snp_nonsynonymous	641399	C	AU5341	E-07	NODE_2
	"""
	t = import_table(string, index = ["seq id", "position"])
	return t


@pytest.fixture
def deletion_group() -> pandas.DataFrame:
	test_deletion_group = pandas.DataFrame(
		{
			'alt':              ['G', 'G'],
			'aminoAlt':         [nan, nan],
			'aminoRef':         [nan, nan],
			'annotation':       [nan, nan],
			'codonAlt':         [nan, nan],
			'codonRef':         [nan, nan],
			'description':      ['large deletion', 'large deletion'],
			'evidence':         ['MC', 'MC'],
			'gene':             ['parA-pQBR0478', 'parA-pQBR0478'],
			'index':            [0, 0],
			'locusTag':         ['[pQBR0001]–[pQBR0478]', '[pQBR0001]–[pQBR0478]'],
			'mutation':         ['D425,094 bp', 'D425,094 bp'],
			'mutationCategory': ['large_deletion', 'large_deletion'],
			'position':         [1, 1],
			'ref':              ['GATCTGCATCTGCAAAGGAGCCAGTCATGTCACCGAAGAAACGCACCCACAAGCCTGCGAAGGTCATTGTCATCGAGAACCAGAAGGGCGGCGTTGGCAA',
				'GATCTGCATCTGCAAAGGAGCCAGTCATGTCACCGAAGAAACGCACCCACAAGCCTGCGAAGGTCATTGTCATCGAGAACCAGAAGGGCGGCGTTGGCAA'],
			'sampleId':         ['SC_58SM', 'S1_58BA'],
			'sampleName':       ['SC_58SM', 'S1_58BA'],
			'seq id':           ['NC_009444', 'NC_009444']
		})
	return test_deletion_group


@pytest.mark.parametrize(
	"values,expected",
	[
		(["12", "13.4", '33'], 58.4 / 3),
		([], 0.0),
		(["", "", ""], 0.0)
	]
)
def test_calculate_average_value(values, expected):
	assert generate_comparison_table._calculate_average_value(values) == expected


def test_parse_mutation_group_snp(snp_group):
	unique_samples = ["S2_58BA", "SC_58SM", "S1_58BA"]
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

	output = generate_comparison_table.parse_mutation_group(snp_group, unique_samples, ref_col, alt_col)
	assert output == expected_output


def test_parse_mutation_group_large_deletion(deletion_group):
	unique_samples = ["S2_58BA", "SC_58SM", "S1_58BA"]
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

	output = generate_comparison_table.parse_mutation_group(deletion_group, unique_samples, ref_col, alt_col)
	assert output == expected_output




def test_extract_string_from_group(mutation_group):
	result = generate_comparison_table._extract_string_from_group(mutation_group, generate_comparison_table.IsolateTableColumns.annotation)
	expected = "L701M (CTT-ATG)"
	assert result == expected


def test_parse_mutation_group(mutation_group):
	expected = {
		'seq id':              'NODE_2',
		'position':            641399,
		'annotation':          "L701M (CTT-ATG)",
		'description':         "amino acid adenylation protein",
		'gene':                'BPFJKKCJ_01788 -',
		'E-01':                "A",
		'E-02':                "A",
		'E-05':                "A",
		'E-06':                "A",
		'E-07':                "A",
		'E-10':                "C",
		'presentIn':           5,
		'presentInAllSamples': False,
		'ref':                 'C',
		'mutationCategory':    'snp_nonsynonymous'
	}
	unique_samples = ['E-01', 'E-02', 'E-05', 'E-06', 'E-07', 'E-10']

	result = generate_comparison_table.parse_mutation_group(mutation_group, unique_samples = unique_samples, ref_col = 'ref', alt_col = 'alt')
	result.pop('locusTag')  # nan cannot be compared against itself.
	assert result == expected


def test_get_relevant_columns():
	assert ('ref', 'alt') == generate_comparison_table._get_relevant_columns('base')

	assert ('aminoRef', 'aminoAlt') == generate_comparison_table._get_relevant_columns('amino')
	assert ('codonRef', 'codonAlt') == generate_comparison_table._get_relevant_columns('codon')


if __name__ == "__main__":
	pass
