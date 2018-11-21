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

def _calculate_average_value(values: List[str])->float:
	values = [i for i in values if i]
	average = sum(float(i) for i in values) / len(values)
	return average

def flatten_mutation_group(unique_samples: List[str], group: pandas.DataFrame) -> Dict[str, Any]:
	ignore = [
		'aaAlt', 'aaRef', 'codonAltSeq', 'codonRefSeq', 'Sample',
		'aaPosition', 'codonNumber', 'coverage', 'polymorphismFrequency'
	]

	row = {k: '|'.join(str(i) for i in group[k].unique() if not pandas.isna(i)) for k in group.columns if k not in (unique_samples + ignore)}

	row['mutationCount'] = len(group)
	row['averageQuality'] = _calculate_average_value(row['quality'].split('|'))
	row['averageDepth'] = _calculate_average_value(row['readDepth'].split('|'))
	row.pop('quality')
	row.pop('readDepth')

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
	#breseq_run_folder = Path("/media/cld100/FA86364B863608A1/Users/cld100/Storage/TravisanoBreseq/")
	breseq_run_folder = Path("/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/lipuma/pipeline_output")

	# df.to_excel(breseq_run_folder / "comparison_table.xlsx")
	sample_map = {'AU0075': 'F-01',
		 'AU0106': 'F-02',
		 'AU0201': 'F-03',
		 'AU0300': 'F-04',
		 'AU0465': 'B-08',
		 'AU10040': 'E-16',
		 'AU1051': 'B-12',
		 'AU1055': 'B-14',
		 'AU1056': 'B-15',
		 'AU1057': 'B-13',
		 'AU1064': 'A-03',
		 'AU10643': 'A-18',
		 'AU11001': 'E-18',
		 'AU11002': 'E-17',
		 'AU11091': 'E-19',
		 'AU1141': 'A-04',
		 'AU1142': 'A-05',
		 'AU1143': 'A-06',
		 'AU11478': 'E-20',
		 'AU12824': 'A-19',
		 'AU13363': 'E-21',
		 'AU13364': 'E-22',
		 'AU14286': 'F-05',
		 'AU15033': 'F-06',
		 'AU15410': 'A-20',
		 'AU15487': 'E-23',
		 'AU15488': 'E-24',
		 'AU1581': 'E-01',
		 'AU1746': 'A-07',
		 'AU1836': 'E-02',
		 'AU18616': 'A-21',
		 'AU20364': 'A-22',
		 'AU20865': 'F-09',
		 'AU20866': 'F-08',
		 'AU21755': 'F-10',
		 'AU23407': 'F-11',
		 'AU23516': 'F-07',
		 'AU2427': 'A-08',
		 'AU2428': 'A-09',
		 'AU25990': 'A-23',
		 'AU25991': 'A-24',
		 'AU2986': 'A-10',
		 'AU31639': 'F-12',
		 'AU33869': 'F-13',
		 'AU3415': 'E-03',
		 'AU3416': 'E-04',
		 'AU35919': 'F-14',
		 'AU36284': 'F-15',
		 'AU36474': 'F-16',
		 'AU36973': 'F-17',
		 'AU3739': 'B-22',
		 'AU3740': 'B-23',
		 'AU3741': 'B-21',
		 'AU37865': 'F-18',
		 'AU3827': 'A-12',
		 'AU3828': 'A-11',
		 'AU4359': 'B-24',
		 'AU4381': 'E-05',
		 'AU4993': 'E-06',
		 'AU5341': 'E-07',
		 'AU6015': 'A-13',
		 'AU6319': 'A-14',
		 'AU6667': 'A-15',
		 'AU6668': 'A-16',
		 'AU6936': 'E-08',
		 'AU7263': 'E-09',
		 'AU7574': 'E-11',
		 'AU7575': 'E-12',
		 'AU7576': 'E-13',
		 'AU7577': 'E-10',
		 'AU8675': 'A-17',
		 'AU8821': 'E-14',
		 'AU9400': 'E-15',
		 'SC1128': 'B-01',
		 'SC1129': 'B-02',
		 'SC1145': 'B-03',
		 'SC1209': 'B-04',
		 'SC1210': 'B-05',
		 'SC1211': 'B-06',
		 'SC1339': 'B-07',
		 'SC1360': 'A-01',
		 'SC1371': 'B-09',
		 'SC1392': 'B-10',
		 'SC1400': 'B-11',
		 'SC1402': 'A-02',
		 'SC1407': 'B-16',
		 'SC1419': 'B-17',
		 'SC1420': 'B-18',
		 'SC1421': 'B-19',
		 'SC1435': 'B-20'}
	variant_df, coverage_df, junction_df = parse_breseqset(breseq_run_folder, whitelist = sample_map.keys())
	comparison_df = generate_snp_comparison_table(variant_df)
	column_whitelist:List[str] = [
		"alt",	"aminoAlt",	"aminoRef",	"annotation",	"description",
		"evidence",	"gene",	"isolates",	"locusTag",	"mutation",
		"mutationCategory",	"mutationCount",	"presentInNSamples",	"quality",	"readDepth",
		"ref",	"variantType",	"presentInAllSamples"
]
	comparison_df.columns = [sample_map.get(i, i) for i in comparison_df.columns]
	sorted_keys = sorted(
		[i for i in comparison_df.columns if '-' in i],
		key = lambda s: (s[0],int(s.split('-')[-1]))
	)
	comparison_df = comparison_df[sorted_keys + column_whitelist]

	tables = {
		'comparison': comparison_df,
		'variant':    variant_df.reset_index(),
		'coverage':   coverage_df.reset_index(),
		'junction':   junction_df.reset_index()
	}

	save_isolate_table(tables, breseq_run_folder / "breseq_table.xlsx")
