from pathlib import Path
from typing import Any, Dict, List, Tuple
from breseqparser.isolate_parser import IsolateTableColumns
import pandas


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


def generate_snp_comparison_table(breseq_table: pandas.DataFrame) -> pandas.DataFrame:
	"""
		Generates a table with sample alt sequences represented by columns
	Parameters
	----------
	breseq_table:TableIO

	Returns
	-------
	pandas.DataFrame
	"""
	if 'filterOut' in breseq_table:
		breseq_table = breseq_table[~breseq_table['filterOut']]
	comparison_table = list()
	groups: List[Tuple[str, pandas.DataFrame]] = breseq_table.groupby(by = [IsolateTableColumns.sequence_id, IsolateTableColumns.position])
	unique_samples = list(breseq_table['sampleId'].unique())

	for key, group in groups:
		number_of_isolates = len(group)
		mutation_data = flatten_mutation_group(unique_samples, group)
		reference_sequence = mutation_data[IsolateTableColumns.ref]
		if reference_sequence == 'nan' or pandas.isna(reference_sequence):
			reference_sequence = ""

		sample_sequences: Dict[str, str] = {k: reference_sequence for k in unique_samples}

		for index, sample in group.iterrows():
			sample_id = sample[IsolateTableColumns.sample_id]
			alt_seq = sample[IsolateTableColumns.alt]

			if pandas.isna(alt_seq) or str(alt_seq) == 'nan':
				alt_seq = reference_sequence

			sample_sequences[sample_id] = alt_seq

		# Check if all sequences are nan, which occurs for large deletions.
		alternate_sequences = set(sample_sequences.values())
		if len(alternate_sequences) == 1 and reference_sequence in alternate_sequences:
			continue

		mutation_group = {**mutation_data, **sample_sequences}
		mutation_group['isolates'] = number_of_isolates
		comparison_table.append(mutation_group)

	df = pandas.DataFrame(comparison_table)

	# Add a filter for variants that appear in all samples.
	total_samples = breseq_table[IsolateTableColumns.sample_id].nunique()
	df['presentInAllSamples'] = df['isolates'] == total_samples
	return df


if __name__ == "__main__":
	from breseqset_parser import parse_breseqset
	from file_generators import save_isolate_table

	breseq_run_folder = Path("/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/lipuma/pipeline_output")

	variant_df, coverage_df, junction_df = parse_breseqset(breseq_run_folder)
	comparison_df = generate_snp_comparison_table(variant_df)

	tables = {
		'comparison': comparison_df,
		'variant':    variant_df.reset_index(),
		'coverage':   coverage_df.reset_index(),
		'junction':   junction_df.reset_index()
	}

	save_isolate_table(tables, breseq_run_folder / "breseq_table.xlsx")
