import pandas
from pathlib import Path
from typing import Union


def calculate_basic_statistics(comparison_table: pandas.DataFrame):
	"""
		Calculates basic statistics related to the mutations found in the comparison table.
	Parameters
	----------
	comparison_table

	Returns
	-------

	"""
	sample_columns = [i for i in comparison_table.columns if '-' in i]
	reference = comparison_table['ref']
	mutation_category = comparison_table['mutationCategory'].copy()

	mutation_category_table = list()
	for sample_column in sample_columns:
		sample = comparison_table[sample_column]
		sample_variants = sample != reference

		sample_mutation_category = mutation_category[sample_variants].value_counts().to_dict()
		sample_mutation_category['sampleName'] = sample_column
		mutation_category_table.append(sample_mutation_category)

	mutation_category_table = pandas.DataFrame(mutation_category_table)
	mutation_category_table = mutation_category_table.fillna(0.0)
	try:
		mutation_category_table['dN/dS'] = mutation_category_table['snp_nonsynonymous'] / mutation_category_table['snp_synonymous']
	except KeyError:
		mutation_category_table['dN/dS'] = 0
	return mutation_category_table


def parse_metadata_table(metadata_table: Union[Path, pandas.DataFrame]) -> pandas.DataFrame:
	if not isinstance(metadata_table, pandas.DataFrame):
		metadata_table = pandas.read_excel(metadata_table)

	metadata_table = metadata_table[["RepositoryNumber", "groupId", "Category", "CultureDate"]]
	metadata_table.columns = ["sampleID", "sampleName", "patient", "date"]
	metadata_table = metadata_table.set_index("sampleName")
	return metadata_table


def random_color() -> str:
	import random
	r = random.randint(0, 255)
	g = random.randint(0, 255)
	b = random.randint(0, 255)
	string = f"#{r:>02X}{g:>02X}{b:>02X}"
	return string


def plot_variants(mutations: pandas.DataFrame, metadata: pandas.DataFrame, output: Path = None):
	colormap = {i: random_color() for i in metadata['patient'].unique()}
	print(mutations['sampleName'].tolist())
	sample_metadata = metadata.loc[mutations['sampleName'].values]
	print(sample_metadata.to_string())
	import matplotlib.pyplot as plt
	dates = sample_metadata['date']
	dates = [pandas.to_datetime(i) for i in dates]
	colormap = [colormap.get(i, "#FF0000") for i in sample_metadata['patient'].values]

	snps = mutations['snp_nonsynonymous']

	fig, ax = plt.subplots(figsize = (10, 5))

	ax.scatter(dates, snps, c = colormap)
	if output:
		plt.savefig(str(output))


if __name__ == "__main__":
	folder = Path("/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/lipuma/sibling_pair_a/")
	metadata_filename = Path("/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/lipuma/docs/merged_table.xlsx")
	table_filename = folder / "SC1360.variants.xlsx"
	table = pandas.read_excel(table_filename)
	output_filename = table_filename.with_suffix(".basicstatistics.tsv")
	metadata = parse_metadata_table(metadata_filename)
	result = calculate_basic_statistics(table)
	plot_variants(result, metadata, table_filename.with_suffix('.png'))
	result.to_csv(output_filename, sep = '\t', index = False)
