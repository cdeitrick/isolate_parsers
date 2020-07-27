from typing import Dict, List, Tuple, Union, Iterable
from loguru import logger
import pandas

from isolateparser.resultparser import IsolateTableColumns


def _calculate_average_value(values: List[str]) -> float:
	""" Calculates the average value of a `|` delimited list of values."""
	values = [i for i in values if i]
	try:
		average = sum(float(i) for i in values) / len(values)
	except ZeroDivisionError:
		average = 0
	return average


def _extract_string_from_group(group: pandas.DataFrame, column: str) -> str:
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

def extract_alt_from_group(group:pandas.DataFrame)->str:

	unique_alts = list(group[IsolateTableColumns.alt].unique())
	# alt may be NaN
	unique_alts = [i for i in unique_alts if isinstance(i, str)]
	alts = "|".join(unique_alts)

	return alts
def parse_mutation_group(group: pandas.DataFrame, unique_samples: List[str], ref_col: str, alt_col: str) -> Dict[str, Union[int, float, str]]:
	# Get a list of all columns that should be identical throughout the mutational group.
	static_columns = [
		IsolateTableColumns.sequence_id, IsolateTableColumns.position,
		IsolateTableColumns.locus_tag,
		IsolateTableColumns.gene, IsolateTableColumns.mutation_category
	]
	# _validate_mutation_group(group, static_columns)
	# Annotation depends on the 'alt' sequence.

	# Retrieve the values for the static columns from the first row.
	first_row = group.reset_index().iloc[0]
	annotation = _extract_string_from_group(group, IsolateTableColumns.annotation)
	description = _extract_string_from_group(group, IsolateTableColumns.description)
	if not annotation: annotation = first_row[IsolateTableColumns.mutation]

	static_data = first_row.reindex(static_columns)
	static_data = static_data.to_dict()
	static_data[IsolateTableColumns.annotation] = annotation
	static_data[IsolateTableColumns.description] = description

	# Large deletions do not have an annotation.
	# Replace with the text in the `mutation` field.
	reference = first_row[ref_col]

	static_data[ref_col] = reference

	static_data[alt_col] = extract_alt_from_group(group)
	for index, row in group.iterrows():
		sample_name = row[IsolateTableColumns.sample_name]
		static_data[sample_name] = row[alt_col]

	# Add the reference value for all samples not part of this mutation group.

	for sample_name in unique_samples:
		if sample_name not in static_data:
			static_data[sample_name] = 0.0 if alt_col == 'frequency' else reference

	static_data['presentInAllSamples'] = len(group) == len(unique_samples)
	static_data['presentIn'] = len(group)
	# When the samples were populations add a column that specifies the alternate value for this group.
	# Use the group rather than the first row in case there are multiple alternate bases
	if alt_col == 'frequency':
		elements = list(group[IsolateTableColumns.alt].unique())
		elements = [(i if isinstance(i, str) else "N/A") for i in elements]
		static_data['alt'] = "|".join(elements)
	# static_data[IsolateTableColumns]

	return static_data


def _get_relevant_columns(by: str, is_population:bool) -> Tuple[str, str]:
	if by == 'base':
		reference_column = IsolateTableColumns.ref
		# Hard coded for now because IsolateTableColumns has the column labels before they were renamed.
		alternate_column = IsolateTableColumns.alt if not is_population else 'frequency'
	elif by == 'amino':
		reference_column = IsolateTableColumns.reference_amino
		alternate_column = IsolateTableColumns.alternate_amino
	elif by == 'codon':
		reference_column = IsolateTableColumns.reference_codon
		alternate_column = IsolateTableColumns.alternate_codon
	else:
		raise ValueError(f"Invalid option: '{by}', expected one of 'base', 'amino', 'codon'")

	return reference_column, alternate_column

def apply_cds_annotations(df:pandas.DataFrame)->pandas.DataFrame:
	""" Genomes downloaded from the ncbi website include a translated_cds file, which can be used to annotated a denovo assembly.
		The resulting annotations (saved in the 'description' field) include a lot of useful metadata in the form of [`key`=`value`].
	"""
	import re

	def _apply(pattern:str, sequence:Iterable)->List[str]:
		result = list()
		for element in sequence:

			match = re.search(pattern, element)

			if match:
				result.append(match.group(1))
			else:
				result.append(element)
		return result

	locus_tag_pattern = "locus_tag=([^\]]+)"
	gene_pattern = "protein=(.+?)\]"


	new_locus_tags = _apply(locus_tag_pattern, df['description'].tolist())
	new_genes = _apply(gene_pattern, df['description'].tolist())

	df['geneOld'] = df['gene'].values
	df['locusTagOld'] = df['locusTag'].values
	df['locusTag'] = new_locus_tags
	df['gene'] = new_genes

	return df

def check_if_population(table:pandas.DataFrame)->bool:
	""" Checks whether the 'frequency' column is non-NaN, indicating that this was a population."""
	freq = table['frequency']

	# If the `frequency` column is entirely NaN, then there will only be one unique value.
	unique = freq.unique()

	return len(unique) != 1

def generate_snp_comparison_table(breseq_table: pandas.DataFrame, by: str,reference_sample: str = None) -> pandas.DataFrame:
	"""
		Generates a table with sample alt sequences represented by columns.
	Parameters
	----------
	breseq_table:pandas.DataFrame
		The concatenated variant tables for all samples.
	by: {'base', 'codon', 'amino'}
		Indicates which reference to use.
	filter_table: bool
		Indicates whether to filter out mutations which fail certain filters.

	reference_sample:str
		Label of the sample to use as a reference. If given, a new column will be added to the table indicating if the reference sample also
		contained the variant.

	Returns
	-------
	pandas.DataFrame
	"""

	unique_samples = list(breseq_table[IsolateTableColumns.sample_name].unique())
	is_population = check_if_population(breseq_table)
	reference_column, alternate_column = _get_relevant_columns(by, is_population)

	_group_by = ['seq id', 'position', 'mutationCategory']
	position_groups: List[Tuple[str, pandas.DataFrame]] = breseq_table.groupby(by = _group_by)

	comparison_table = list()
	for key, group in position_groups:
		result = parse_mutation_group(group, unique_samples, reference_column, alternate_column)
		comparison_table.append(result)
	df = pandas.DataFrame(comparison_table)
	if df.empty:
		message = f"The comparison table could not be created, and is empty. This may be due to the `mutationCategory` column being empty."
		logger.warning(message)

	# Add a column indicating if the reference sample contained the variant. This only applies if the reference sample is known.
	if reference_sample and reference_sample in df.columns:
		df['inReference'] = ((df[reference_sample] != df[reference_column]) & (df[IsolateTableColumns.mutation_category] != "large_deletion"))
	# Check if the description field came from a translated cds file.

	#df = apply_cds_annotations(df)


	return df
