from pathlib import Path
from typing import Optional, Union, List
import io
import pandas


def _import_table_from_path(filename: Path, sheet_name: Optional[str] = None, index: Optional[str] = None) -> pandas.DataFrame:
	""" Imports a file as a pandas.DataFrame. Infers filetype from the filename extension/suffix.
	"""
	if filename.suffix in {'.xls', '.xlsx'}:
		data: pandas.DataFrame = pandas.read_excel(str(filename), sheet_name = sheet_name)
	else:
		sep = '\t' if filename.suffix in {'.tsv', '.tab'} else ','
		data: pandas.DataFrame = pandas.read_table(str(filename), sep = sep)

	if index and index in data.columns:
		data = data.set_index(index)

	return data


def _import_table_from_string(string: str, delimiter: Optional[str] = None, index: Optional[str] = None) -> pandas.DataFrame:
	""" Imports a table represented as a basic string object."""
	# Remove unwanted whitespace.
	string = '\n'.join(i.strip() for i in string.split('\n') if i)
	if not delimiter:
		delimiter = '\t' if '\t' in string else ','
	result = pandas.read_table(io.StringIO(string), sep = delimiter, index_col = False)
	if index:
		# Using `index_col` in `read_table()` doesn't work for some reason.
		#result[index] = result[index].astype(str)
		result.set_index(index, inplace = True)
	return result


def import_table(input_table: Union[str, Path], sheet_name: Optional[str] = None, index:Optional[Union[str, List[str]]] = None) -> pandas.DataFrame:
	if isinstance(input_table, Path):
		data = _import_table_from_path(input_table, sheet_name, index)
	else:
		data = _import_table_from_string(input_table, index = index)
	#data = data[sorted(data.columns, key = lambda s: str(s))]
	return data
