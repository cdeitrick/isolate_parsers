""" Generates an aligned fasta file with the concatenated SNPs from the variant table."""
from functools import partial
from pathlib import Path
from typing import Tuple, Optional

import pandas

from breseqparser.isolate_parser import IsolateTableColumns

SEQUENCE_ID_COLUMN = IsolateTableColumns.sequence_id
POSITION_COLUMN = IsolateTableColumns.position


def generate_reference_sequence(snp_table: pandas.DataFrame, reference_column: str) -> pandas.Series:
	# indexed_table = snp_table.reset_index().set_index(keys = [SEQUENCE_ID_COLUMN, POSITION_COLUMN], verify_integrity = False)
	unindexed_table = snp_table.reset_index()
	groups = unindexed_table.groupby(by = [SEQUENCE_ID_COLUMN, POSITION_COLUMN])
	reference_data = list()

	for (seq_id, position), group in groups:
		_unique_values = group[reference_column].unique()
		if len(_unique_values) != 1:
			print(f"Found {_unique_values} values.")
			print(group.to_string())
			raise ValueError
		reference_sequence = group[reference_column].iloc[0]
		row = {
			SEQUENCE_ID_COLUMN: seq_id,
			POSITION_COLUMN:    position,
			'reference':        reference_sequence
		}
		reference_data.append(row)

	reference_table = pandas.DataFrame(reference_data).set_index([SEQUENCE_ID_COLUMN, POSITION_COLUMN])
	reference_table = reference_table.sort_index()
	return reference_table['reference']


def _parse_sample_group(sample_name: str, group: pandas.DataFrame, reference_sequence: pandas.Series, alt_column: str) -> pandas.Series:
	""" Converts a group of rows from the output of `parse_breseqset` corresponding to a single sample.
		Parameters
		----------
		sample_name: str
			The Id for this sample.
		reference_sequence:pandas.Series
			A pandas.Series object with all reference positions available for the full breseq set.
			- Index -> (`seq id`, `position`)
			- values _> str
		group: pandas.DataFrame
			A subset of the breseqset table corresponding to a single sample.
		Returns
		-------
		pandas.Series
			- Index -> `seq id`, `position`
			- Values-> `alt`
			- name  -> `sample_name`
	"""
	# Use `seq id` and `position` as indicies. These should form a unique tuple that can be mapped back to the reference.
	group = group.set_index([SEQUENCE_ID_COLUMN, POSITION_COLUMN])

	# We are only interested in the `alt` sequences.
	sample_alt: pandas.Series = group[alt_column]
	sample_alt.name = sample_name

	# Align the sample sequence to the reference.
	sample_alt, _ = sample_alt.align(reference_sequence, join = 'right')

	# The sample group only contains positions that differ from the reference. Should align to the full reference (ref for all sequences)
	sample_alt = sample_alt.where(sample_alt.notna(), other = reference_sequence)
	return sample_alt


def _convert_combined_table_to_aligned_table(snp_table: pandas.DataFrame, reference_sequence: pandas.Series, alt_col: str, reference_label:str) -> pandas.DataFrame:
	"""
		Converts the combined breseq table generated by the breseqset parser into a `pandas.DataFrame` object.
	Parameters
	----------
	snp_table: pandas.DataFrame
		The snp table generated by the breseqset parser.

	Returns
	-------
	pandas.DataFrame
		A table with rows corresponding to a single sample and columns corresponding to (`seq id`, `position`) indicies.
	"""
	partial_parse_sample = partial(
		_parse_sample_group,
		reference_sequence = reference_sequence,
		alt_column = alt_col
	)
	groups = snp_table.groupby(by = IsolateTableColumns.sample_name)

	sample_alts = [reference_sequence] + [partial_parse_sample(name, group) for name, group in groups]
	df: pandas.DataFrame = pandas.concat(sample_alts, axis = 1)

	# Filter out variants that are present in the reference sample
	if reference_label:
		reference_sample_variants = df['reference'] != df[reference_label]
		df = df[~reference_sample_variants]

	# It is easier to iterate over rows rather than columns, so transpose the dataframe such that rows correspond to samples.
	df = df.transpose()
	return df


def _validate_variant_table(variant_table: pandas.DataFrame, by: str, ref_col: str, alt_col: str) -> None:
	""" Raises an error if the variant table has errors."""
	# Since this only concerns SNPs, make sure the alt and ref sequences are single characters
	unique_reference_bases = set(variant_table[ref_col].unique())
	unique_alternate_bases = set(variant_table[alt_col].unique())
	if by == 'base':
		try:
			assert unique_reference_bases <= {'A', 'C', 'T', 'G', 'N'}
			assert unique_alternate_bases <= {'A', 'C', 'T', 'G', 'N'}
		except AssertionError:
			print("Found invalid characters!")
			print("reference", unique_reference_bases)
			print("alternate", unique_alternate_bases)
			raise ValueError


def _get_relevant_columns(by: str) -> Tuple[str, str]:
	""" Retrieves the relevant reference and alternate columns for the given fasta category."""
	if by == 'codon':
		reference_column = IsolateTableColumns.reference_codon
		alternate_column = IsolateTableColumns.alternate_codon
	elif by == 'amino':
		reference_column = IsolateTableColumns.reference_amino
		alternate_column = IsolateTableColumns.alternate_amino
	elif by == 'base':
		reference_column = IsolateTableColumns.ref
		alternate_column = IsolateTableColumns.alt
	else:
		message = f"Invalid fasta option '{by}', expected one of 'codon', 'amino', 'base'"
		raise ValueError(message)

	return reference_column, alternate_column


def write_fasta_file(df: pandas.DataFrame, filename: Path) -> Path:
	""" Writes each row of a dataframe to a fasta file. Each row should correspond to a specific sample and should be indexed by sample name or id."""
	with filename.open('w') as fasta_file:
		for index, row in df.iterrows():
			seq = "".join(row.tolist())
			fasta_file.write(f">{index}\n{seq}\n")
	return filename

def _filter_out_variants_in_reference(aligned_table:pandas.DataFrame, reference_label:str):
	pass

def generate_fasta_file(variant_table: pandas.DataFrame, filename: Path, by: str = 'codon', reference_label:Optional[str] = None) -> pandas.DataFrame:
	"""Converts the variant table generated from the breseqset parser into a fasta file."""
	reference_column, alternate_column = _get_relevant_columns(by)

	# Make sure `filename` is a Path
	filename = Path(filename)
	table = variant_table.reset_index()

	# We only care about snps.
	accepted_mutation_categories = ['snp_nonsynonymous']#, 'snp_synonymous']
	table = table[table[IsolateTableColumns.mutation_category].isin(accepted_mutation_categories)]
	groups = table.groupby(by = [IsolateTableColumns.sample_name, IsolateTableColumns.mutation_category])
	samples = dict()
	for (isolate_name, mutation_category), group in groups:
		if isolate_name in samples:
			samples[isolate_name][mutation_category] = len(group)
		else:
			samples[isolate_name] = {mutation_category:len(group)}

	_validate_variant_table(table, by, reference_column, alternate_column)

	reference_sequence: pandas.Series = generate_reference_sequence(table, reference_column)
	df = _convert_combined_table_to_aligned_table(table, reference_sequence, alternate_column, reference_label)

	# Write the dataframe to a fasta file.
	write_fasta_file(df, filename)

	return df


if __name__ == "__main__":
	from breseqparser.isolate_parser import parse_breseq_isolate

	_filename = Path("/home/cld100/Documents/github/isolate_parsers/tests/data/Clonal_Output/")
	_snp_table, *_ = parse_breseq_isolate(_filename, 'testIsolate', 'testIsolate')

	print(_snp_table.to_string())
	generate_fasta_file(_snp_table, _filename.parent / "test_fasta.codon.fasta", by = 'codon')
	generate_fasta_file(_snp_table, _filename.parent / "test_fasta.base.fasta", by = 'base')
