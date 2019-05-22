
from typing import Dict, List, Tuple, Union
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




def _extract_string_from_group(group:pandas.DataFrame, column:str)->str:
	"""
		Essentially concatenates the annotations for each of the samples in the group. This will typically be the same for all variants,
		but is assumed otherwise just in case.
	Parameters
	----------
	group: pandas.DataFrame

	Returns
	-------

	"""
	annotation_values = group[column].tolist()
	# Remove duplicate and missing annotations. DataFrames save missing values as math.nan
	annotation_set = {i for i in annotation_values if isinstance(i, str)}
	# Combine and merge into a string
	annotation_string = "|".join(annotation_set)
	return annotation_string

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
	annotation = _extract_string_from_group(group, IsolateTableColumns.annotation)
	description = _extract_string_from_group(group, IsolateTableColumns.description)
	if not annotation: annotation = first_row[IsolateTableColumns.mutation]
	static_data = first_row[static_columns].to_dict()
	static_data[IsolateTableColumns.annotation] = annotation
	static_data[IsolateTableColumns.description] = description

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

def generate_snp_comparison_table(breseq_table: pandas.DataFrame, by: str, filter_table:bool = False, reference_sample:str = None) -> pandas.DataFrame:
	"""
		Generates a table with sample alt sequences represented by columns.
	Parameters
	----------
	breseq_table:pandas.DataFrame
		The concatenated variant tables for all samples.
	by: {'base', 'codon', 'amino'}
		Indicates which reference to use.
	reference_sample:str
		Label of the sample to use as a reference. If given, a new column will be added to the table indicating if the reference sample also
		contained the variant.

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

	# Add a column indicating if the reference sample contained the variant. This only applies if the reference sample is known.
	if reference_sample and reference_sample in df.columns:
		df['inReference'] = df[reference_sample] != df[reference_column]
	return df

