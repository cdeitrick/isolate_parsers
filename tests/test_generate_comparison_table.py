import pandas
import pytest

import dataio
from isolateparser.generate import generate_comparison_table


@pytest.fixture
def mutation_group_clones() -> pandas.DataFrame:
	string = """
		alt	aminoAlt	aminoRef	annotation	codonAlt	codonRef	description	evidence	gene	locusTag	mutation	mutationCategory	position	ref	sampleId	sampleName	seq id
		A	M	L	L701M (CTT-ATG)	ATG	CTT	amino acid adenylation protein	RA	BPFJKKCJ_01788 -		C-A	snp_nonsynonymous	641399	C	AU1581	E-01	NODE_2
		A	M	L	L701M (CTT-ATG)	ATG	CTT	amino acid adenylation protein	RA	BPFJKKCJ_01788 -		C-A	snp_nonsynonymous	641399	C	AU1836	E-02	NODE_2
		A	M	L	L701M (CTT-ATG)	ATG	CTT	amino acid adenylation protein	RA	BPFJKKCJ_01788 -		C-A	snp_nonsynonymous	641399	C	AU4381	E-05	NODE_2
		A	M	L	L701M (CTT-ATG)	ATG	CTT	amino acid adenylation protein	RA	BPFJKKCJ_01788 -		C-A	snp_nonsynonymous	641399	C	AU4993	E-06	NODE_2
		A	M	L	L701M (CTT-ATG)	ATG	CTT	amino acid adenylation protein	RA	BPFJKKCJ_01788 -		C-A	snp_nonsynonymous	641399	C	AU5341	E-07	NODE_2
	"""
	t = dataio.import_table(string, index = ["seq id", "position"])
	return t


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


def test_extract_string_from_group(mutation_group_clones):
	result = generate_comparison_table._extract_string_from_group(mutation_group_clones, generate_comparison_table.IsolateTableColumns.annotation)
	expected = "L701M (CTT-ATG)"
	assert result == expected


def test_get_relevant_columns_for_clones():
	assert ('ref', 'alt') == generate_comparison_table._get_relevant_columns('base', is_population = False)

	assert ('aminoRef', 'aminoAlt') == generate_comparison_table._get_relevant_columns('amino', is_population = False)
	assert ('codonRef', 'codonAlt') == generate_comparison_table._get_relevant_columns('codon', is_population = False)


def test_get_relevant_colums_for_populations():
	assert ('ref', 'frequency') == generate_comparison_table._get_relevant_columns('base', True)


if __name__ == "__main__":
	pass
