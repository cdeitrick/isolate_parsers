from pathlib import Path
from typing import Dict

import pandas


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
		df.to_excel(writer, sheet_label, index = include_index)
	return filename
