import pandas

from breseqset_parser import parse_breseqset, save_isolate_table


def generate_comparison_table(variant_df: pandas.DataFrame) -> pandas.DataFrame:
	"""

	Parameters
	----------
	variant_df

	Returns
	-------
	"""
	variant_df = variant_df.reset_index()
	variant_df = variant_df.set_index(['position'])
	variant_df.pop('quality')
	variant_df.pop('readDepth')
	# print(variant_df.to_string())

	groups = variant_df.groupby(by = 'Sample')
	sample_df = variant_df.pivot(columns = 'Sample', values = 'alt')
	print(sample_df.to_string())
	# merged_df = variant_df.merge(sample_df, left_index = True, right_index = True, how = 'left')
	merged_df = variant_df.join(sample_df)
	print(merged_df.to_string())


from pathlib import Path
import pandas
from typing import TypeVar, List, Tuple, Any, Dict

TableIO = TypeVar('IO', Path, pandas.DataFrame)
GroupType = List[Tuple[str, pandas.DataFrame]]


def flatten_mutation_group(unique_samples: List[str], group: pandas.DataFrame) -> Dict[str, Any]:
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

	"""
	if 'filterOut' in breseq_table:
		breseq_table = breseq_table[~breseq_table['filterOut']]
	comparison_table = list()
	groups: GroupType = breseq_table.groupby(by = ['seq id', 'position'])
	unique_samples = list(breseq_table['Sample'].unique())

	for key, group in groups:
		number_of_isolates = len(group)
		mutation_data = flatten_mutation_group(unique_samples, group)
		reference_sequence = mutation_data['ref']
		if reference_sequence == 'nan' or pandas.isna(reference_sequence):
			reference_sequence = ""

		sample_sequences: Dict[str, str] = {k: reference_sequence for k in unique_samples}

		for index, sample in group.iterrows():
			sample_id = sample['Sample']
			alt_seq = sample['alt']

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
	total_samples = breseq_table['Sample'].nunique()
	df['presentInAllSamples'] = df['isolates'] == total_samples
	return df


if __name__ == "__main__":
	# breseq_run_folder = Path(__file__).parent.parent /"data"/ "breseq_run"
	breseq_run_folder = Path("/media/cld100/FA86364B863608A1/Users/cld100/Storage/TravisanoBreseq/")
	variant_df, coverage_df, junction_df = parse_breseqset(breseq_run_folder)
	comparison_df = generate_snp_comparison_table(variant_df)
	# df.to_excel(breseq_run_folder / "comparison_table.xlsx")
	tables = {
		'comparison': comparison_df,
		'variant':    variant_df.reset_index(),
		'coverage':   coverage_df.reset_index(),
		'junction':   junction_df.reset_index()
	}

	save_isolate_table(tables, breseq_run_folder / "breseq_table.xlsx")
