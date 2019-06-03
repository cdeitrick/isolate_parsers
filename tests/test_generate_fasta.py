from pathlib import Path

import pandas
import pytest

import dataio
# Need to also import private functions
from isolateparser.file_generators import generate_fasta

data_folder = Path(__file__).parent / 'data' / 'Clonal_Output' / 'breseq_output'
gd_filename = data_folder / 'output' / 'evidence' / 'annotated.gd'


@pytest.fixture
def variant_table():
	test_table = [
		{'ref': 'G', 'alt': 'A', 'sampleName': 'S2_58BA', 'seq id': 'NC_012660', 'position': 3959631},
		{'ref': 'A', 'alt': 'C', 'sampleName': 'SC_58SM', 'seq id': 'NC_012660', 'position': 4173231},
		{'ref': 'G', 'alt': 'A', 'sampleName': 'S1_58BA', 'seq id': 'NC_012660', 'position': 3959631}
	]
	return pandas.DataFrame(test_table)


@pytest.fixture
def small_table() -> pandas.DataFrame:
	string = """
	E-21	E-22	E-23	E-24	ref
	C	C	C	C	G
	G	G	G	G	G
	T	T	T	G	T
	G	G	GG	G	G
	A	A	A	A	C
	T	T	T	T	G
	"""
	table = dataio.import_table(string)
	return table


def test_write_fasta_file(variant_table, tmp_path):
	expected = """
	>ref
	GGTGCG
	>E-21
	CGTGAT
	>E-22
	CGTGAT
	>E-23
	CGTGGAT
	>E-24
	CGGGAT
	"""
	expected = "\n".join(j for j in [i.strip() for i in expected.split('\n')] if j) + '\n'
	df = variant_table.transpose()
	df = df.loc[['ref', 'E-21', 'E-22', 'E-23', 'E-24']]  # To make sure the file is ordered correctly.
	temporary_file = tmp_path / "temp.fasta"
	generate_fasta.write_fasta_file(df, temporary_file)
	assert temporary_file.read_text() == expected


def test_filter_variants_in_sample(variant_table):
	expected = """
		E-21	E-22	E-23	E-24	ref
		G	G	G	G	G
		G	G	GG	G	G
	"""
	expected = dataio.import_table(expected)

	result = generate_fasta._filter_variants_in_sample(variant_table, 'E-24', 'ref')
	result = result.reset_index()
	result.pop('index')

	pandas.testing.assert_frame_equal(expected, result)


def test_get_relevant_columns():
	ref_base, alt_base = generate_fasta._get_relevant_columns('base')
	assert 'ref' == ref_base
	assert 'alt' == alt_base

	ref_codon, alt_codon = generate_fasta._get_relevant_columns('codon')
	assert ref_codon == 'codonRef'
	assert alt_codon == 'codonAlt'

	ref_amino, alt_amino = generate_fasta._get_relevant_columns('amino')
	assert ref_amino == 'aminoRef'
	assert alt_amino == 'aminoAlt'


def test_generate_reference_sequence(variant_table):
	test_df = pandas.DataFrame(variant_table)
	truth_reference: pandas.Series = pandas.DataFrame(
		[
			{'seq id': 'NC_012660', 'position': 3959631, 'reference': 'G'},
			{'seq id': 'NC_012660', 'position': 4173231, 'reference': 'A'}
		]
	).set_index(keys = ['seq id', 'position'])['reference']

	test_reference = generate_fasta.generate_reference_sequence(test_df, 'ref')

	pandas.testing.assert_series_equal(truth_reference, test_reference)


def test_parse_sample_group():
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

	output = generate_fasta._parse_sample_group(sample_group[0], sample_group[1], truth_reference, 'alt')

	pandas.testing.assert_series_equal(truth_output, output)
