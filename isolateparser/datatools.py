import pandas


def variant_count(table: pandas.DataFrame, reference: pandas.Series) -> pandas.Series:
	"""
	Counts the number of samples a specific variant occurs in.
	Parameters
	----------
	table: pandas.DataFrame
		The table should be indexed by sequence Id and position, and have samples as columns and varaints as rows.
	reference:pandas.Series
		The reference series used to determine if a given position counts as a variant.
	Returns
	-------
	pandas.Series
		Maps indexed positions to a count of the number of samples a variant occurs in.
	"""

	present_df: pandas.DataFrame = table.apply(lambda s: s != reference)
	return present_df.sum(axis = 1)


def filter_variants_in_all_samples(df: pandas.DataFrame, reference_label: str) -> pandas.DataFrame:
	"""
	Filters out variant which appear in all samples.
	Parameters
	----------
	df: pandas.DataFrame
		A Dataframe where columns correspond to samples and rows correspond to unique variants.
		The table should be indexed by sequence Id and position.
	reference_label: str
		Used to annotate mutations that appear in the reference sample.
	"""

	reference = df.pop(reference_label)
	variants = variant_count(df, reference)
	present_in_all = variants == len(df.columns)
	return df[~present_in_all]
