from pathlib import Path
from typing import Any, Dict, List, Tuple, Union
import math
import pandas

from breseqparser.isolate_parser import IsolateTableColumns


def _calculate_average_value(values: List[str]) -> float:
	""" Calculates the average value of a `|` delimited list of values."""
	values = [i for i in values if i]
	try:
		average = sum(float(i) for i in values) / len(values)
	except ZeroDivisionError:
		average = 0
	return average


def flatten_mutation_group(unique_samples: List[str], group: pandas.DataFrame) -> Dict[str, Any]:
	""" Averages the values of `readDepth` and `readQuality` from all samples."""
	ignore = [
		'aaAlt', 'aaRef', 'codonAltSeq', 'codonRefSeq', 'Sample',
		'aaPosition', 'codonNumber', 'coverage', 'polymorphismFrequency'
	]

	row = {k: '|'.join(str(i) for i in group[k].unique() if not pandas.isna(i)) for k in group.columns if k not in (unique_samples + ignore)}

	row['mutationCount'] = len(group)

	return row


def _validate_mutation_group(group: pandas.DataFrame, static_columns: List[str]) -> pandas.DataFrame:
	"""Determines whther a column contains more than one unique value."""
	for column in static_columns:
		unique_values = group[column].unique()
		if len(unique_values) != 1:
			message = f"Found {unique_values}"
			raise ValueError(message)
	return group


def parse_mutation_group(group: pandas.DataFrame, unique_samples: List[str], ref_col: str, alt_col: str) -> Dict[str, Union[int, float, str]]:
	# Get a list of all columns that should be identical throughout the mutational group.
	static_columns = [
		IsolateTableColumns.sequence_id, IsolateTableColumns.position,
		IsolateTableColumns.locus_tag,
		IsolateTableColumns.gene, IsolateTableColumns.mutation_category
	]
	#_validate_mutation_group(group, static_columns)
	# Annotation depends on the 'alt' sequence.

	# Retrieve the values for the static columns from the first row.
	first_row = group.reset_index().iloc[0]
	annotation = "|".join([i for i in sorted(set(group[IsolateTableColumns.annotation].tolist())) if isinstance(i, str)])
	if not annotation: annotation = first_row[IsolateTableColumns.mutation]
	static_data = first_row[static_columns].to_dict()
	static_data[IsolateTableColumns.annotation] = annotation

	# Large deletions do not have an annotation.
	# Replace with the text in the `mutation` field.
	reference = first_row[ref_col]
	static_data[ref_col] = reference
	_number_of_alternate_samples = len(group)

	for index, row in group.iterrows():
		sample_name = row[IsolateTableColumns.sample_name]
		alternate = row[alt_col]

		static_data[sample_name] = alternate

	for sample_name in unique_samples:
		if sample_name not in static_data:
			static_data[sample_name] = reference

	static_data['presentInAllSamples'] = len(group) == len(unique_samples)
	static_data['presentIn'] = len(group)
	#static_data[IsolateTableColumns]

	return static_data


def _get_relevant_columns(by: str) -> Tuple[str, str]:
	if by == 'base':
		reference_column = IsolateTableColumns.ref
		alternate_column = IsolateTableColumns.alt
	elif by == 'amino':
		reference_column = IsolateTableColumns.reference_amino
		alternate_column = IsolateTableColumns.alternate_amino
	elif by == 'codon':
		reference_column = IsolateTableColumns.reference_codon
		alternate_column = IsolateTableColumns.alternate_codon
	else:
		raise ValueError(f"Invalid option: '{by}', expected one of 'base', 'amino', 'codon'")

	return reference_column, alternate_column


def generate_snp_comparison_table(breseq_table: pandas.DataFrame, by: str, filter_table:bool = False) -> pandas.DataFrame:
	"""
		Generates a table with sample alt sequences represented by columns.
	Parameters
	----------
	breseq_table:pandas.DataFrame
		The concatenated variant tables for all samples.
	by: {'base', 'codon', 'amino'}
		Indicates which reference to use.

	Returns
	-------
	pandas.DataFrame
	"""
	unique_samples = list(breseq_table[IsolateTableColumns.sample_name].unique())
	reference_column, alternate_column = _get_relevant_columns(by)
	if filter_table:
		if 'filterOut' in breseq_table:
			breseq_table = breseq_table[~breseq_table['filterOut']]
	_group_by = [IsolateTableColumns.sequence_id, IsolateTableColumns.position, IsolateTableColumns.mutation_category]
	position_groups: List[Tuple[str, pandas.DataFrame]] = breseq_table.groupby(by = _group_by)

	comparison_table = [
		parse_mutation_group(g, unique_samples, reference_column, alternate_column) for k, g in position_groups
	]
	df = pandas.DataFrame(comparison_table)
	return df


if __name__ == "__main__":
	from breseqset_parser import parse_breseqset

	breseq_run_folder = Path(__file__).parent.parent / "tests" / "data" / "set_output"

	variant_df, coverage_df, junction_df = parse_breseqset(breseq_run_folder)
	comparison_df = generate_snp_comparison_table(variant_df, 'base')

	tables = {
		'comparison': comparison_df,
		'variant':    variant_df.reset_index(),
		'coverage':   coverage_df.reset_index(),
		'junction':   junction_df.reset_index()
	}

	print(comparison_df.to_dict(orient = 'records'))
