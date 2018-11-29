import pandas
from pathlib import Path
from typing import Dict
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
	"""
	writer.close()

	wb = load_workbook(filename)
	ws = wb['junctions']

	merge_columns = [
		c
		for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
		if ws[c + '1'].value in ['Sample', 0, '0', 'freq', 'product', 'score']
	]

	for x in range(0, len(self.junction_table), 2):
		for column in merge_columns:
			cells_to_merge = '{0}{1}:{0}{2}'.format(column, x + 2, x + 3)
			ws.merge_cells(cells_to_merge)

	wb.save(filename)
	"""
	return filename