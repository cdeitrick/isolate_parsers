from pathlib import Path
from typing import Dict

import pandas

from isolateparser.resultparser import IsolateTableColumns

try:
	import openpyxl
except ModuleNotFoundError:
	openpyxl = None


def save_isolate_table(tables: Dict[str, pandas.DataFrame], filename: Path) -> Path:
	"""
		Saves the parsed table as an Excel spreadsheet.
	Parameters
	----------
	tables: Dict[str,pandas.DataFrame]
		A mapping of sheet names to dataframes.
	filename: str, pathlib.Path
		The output file.

	Returns
	-------

	"""
	writer = pandas.ExcelWriter(str(filename))
	include_index = False
	# python 3.5 or 3.6 made all dicts ordered by default, so the sheets will be ordered in the same order they were defined in `tables`
	for sheet_label, df in tables.items():
		if df is None: continue
		df.to_excel(writer, sheet_label, index = include_index)
	writer.save()  # otherwise color_table_cells will not be able to load the file
	# Color in the spreadsheet cells based on whether the sequence differs from the reference.
	if openpyxl is not None:
		try:
			color_table_cells(filename)
		except:
			pass
	return filename


def _get_relevant_columns(worksheet: openpyxl.Workbook):
	# Get the first row
	first_row = worksheet[1]
	# Find the columns that correspond to samples
	excluded_columns = list(IsolateTableColumns) + ['presentIn', 'presentInAllSamples']
	sample_columns = [i.column for i in first_row if i.value not in excluded_columns]
	reference_sequence_column = [i.column for i in first_row if i.value == IsolateTableColumns.ref][0]
	return reference_sequence_column, sample_columns


def color_table_cells(filename: Path):
	""" Colors in the cells of the comparison table to highlight differences between samples.
		Parameters
		----------
		filename: Path
			Path to the excel file. The sheet containing the comparison table should be named 'variant comparison'.
	"""
	workbook = openpyxl.load_workbook(filename = str(filename))

	worksheet = workbook['variant comparison']

	# There is an issue with libreoffice when 'bgColor' is used instead of 'fgColor' where cells are rendered with a black background.
	variant_pattern = openpyxl.styles.PatternFill(fgColor = "FC8D62", fill_type = "solid")
	reference_column_label, sample_column_labels = _get_relevant_columns(worksheet)

	for sample_column_label in sample_column_labels:
		sample_column = worksheet[sample_column_label]
		for cell in sample_column:
			cell_row = cell.row
			if cell_row == 1: continue  # Is the header column
			reference_cell = worksheet[f"{reference_column_label}{cell_row}"]
			is_variant = cell.value != reference_cell.value
			if is_variant:
				cell.fill = variant_pattern
	workbook.save(str(filename.with_suffix('.edited.xlsx')))


if __name__ == "__main__":
	path = Path.home() / "Documents" / "breseq_table.xlsx"
	color_table_cells(path)
