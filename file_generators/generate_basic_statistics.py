import pandas
from file_generators.generate_comparison_table import IsolateTableColumns
def calculate_basic_statistics(comparison_table:pandas.DataFrame):
	"""
		Calculates basic statistics related to the mutations found in the comparison table.
	Parameters
	----------
	comparison_table

	Returns
	-------

	"""
	sample_columns = [i for i in comparison_table.columns if i not in IsolateTableColumns]
	mutation_groups = comparison_table[IsolateTableColumns.mutation_category]
	sample_info = dict()
	for mutation_category, group in mutation_groups:
		total = len(group)
		print(total, mutation_category)



