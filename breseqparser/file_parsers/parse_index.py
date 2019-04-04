from collections import OrderedDict
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union, Optional

import pandas
from bs4 import BeautifulSoup
from unidecode import unidecode

TableType = List[Dict[str, Any]]
DFType = pandas.DataFrame

VariantTableColumnMap = {
	'Sample': 'sample',
	'annotation': 'annotation',
	'description': 'description',
	'freq': 'frequency',
	'gene': 'gene',
	'mutation': 'mutation',
	'position': 'position',
	'seq iq': 'seq id'
}

def extract_sample_name(filename: Path, sample_name:Optional[str] = None):
	# Attempts to infer the sample name from the filename of the index file.
	if sample_name: return sample_name



def to_integer(string: str) -> int:
	""" Converts a string to a number"""
	try:
		string =int(float(string.replace(',', '')))
	except ValueError:
		pass
	return string


def convert_to_dataframe(table: TableType) -> pandas.DataFrame:
	""" Some of the column names may be encoded differently than expected. This method fixes them."""
	table = pandas.DataFrame(table)
	if 'seq\xa0id' in table.columns:
		table['seq id'] = table['seq\xa0id']
		del table['seq\xa0id']
	return table


def _load_index_file(filename: Path) -> BeautifulSoup:
	"""Loads the index file."""
	filename = Path(filename)
	# Use read_bytes to avoid encoding errors.
	soup = BeautifulSoup(filename.read_bytes(), 'lxml')
	return soup


def _extract_variant_table_headers(text: str) -> List[str]:
	begin_snp_header_string = r'<th>evidence</th>'
	end_snp_header_string = '<!-- Item Lines -->'

	begin_snp_header = text.find(begin_snp_header_string)
	end_snp_header = text.find(end_snp_header_string)

	snp_header_full = text[begin_snp_header:end_snp_header]
	# snp_header = snp_header_full[end_snp_header:]

	snp_header_soup = BeautifulSoup(snp_header_full, 'lxml')
	snp_header_soup = [i.text for i in snp_header_soup.find_all('th')]
	return snp_header_soup


def _extract_coverage_and_junction_tables(alph_soup: str) -> Tuple[BeautifulSoup, BeautifulSoup]:
	"""
		Extracts the headers for the snp table, the junction table, and the coverage table.
	Parameters
	----------
	alph_soup: BeautifulSoup
		Equivilent to BeautifulSoup(index.html)

	Returns
	-------
		snp_header, coverage_table, junction_table

	"""

	begin_umc = alph_soup.find(
		'<tr><th align="left" class="missing_coverage_header_row" colspan="11">Unassigned missing coverage evidence</th></tr>')
	end_umc = alph_soup.find(
		'<th align="left" class="new_junction_header_row" colspan="12">Unassigned new junction evidence</th>')
	coverage_string = alph_soup[begin_umc:end_umc]
	junction_string = alph_soup[end_umc:]
	coverage_soup = BeautifulSoup(coverage_string, 'lxml')
	junction_soup = BeautifulSoup(junction_string, 'lxml')
	return coverage_soup, junction_soup


def _extract_index_tables(soup: BeautifulSoup) -> Tuple[
	List[str], BeautifulSoup, BeautifulSoup, BeautifulSoup]:
	"""
		Extracts the relevant tables from the index table.
	Parameters
	----------
	soup: BeautifulSoup
		A BeautifulSoup-parsed version of the index file.

	Returns
	-------
		snp_header, snp_table, coverage_soup, junction_soup
	"""
	soup_text = str(soup)
	snp_header_soup = _extract_variant_table_headers(soup_text)

	normal_table = soup.find_all(attrs = {'class': 'normal_table_row'})
	poly_table = soup.find_all(attrs = {'class': 'polymorphism_table_row'})
	snp_table = normal_table + poly_table

	coverage_soup, junction_soup = _extract_coverage_and_junction_tables(soup_text)

	return snp_header_soup, snp_table, coverage_soup, junction_soup

def _parse_html_row(row:Dict[str,str])->Dict[str,str]:
	"""
		Converts an individual row in one of the html tables to a dictionary.
	Parameters
	----------
	row: Dict[str,str]
		The extracted row from the varaint table present in the index.html file.

	Returns
	-------

	"""
	try:
		#
		row['position'] = to_integer(row['position'])
	except KeyError:
		row['position'] = None
	try:
		row['freq %'] = float(row['freq'][:-1])
	# row.pop('freq')
	except KeyError:
		pass
	# `description` is the column name in the index file.
	if 'javascript' in row['description'] or 'Javascript' in row['description']:
		row['description'] = 'large deletion'

	return row


def _parse_snp_table(sample_name: str, headers: List[str], rows: BeautifulSoup) -> TableType:
	"""
		Parses the SNP table.
	Parameters
	----------
	sample_name: str
		The name of the sample. Usually extracted from the name of the analysis folder.
	headers: List[str]
		Column names for the snp table.
	rows: List
		The rows of the snp table.

	Returns
	-------

	"""
	converted_table = list()
	for tag in rows:
		values = [v.text for v in tag.find_all('td')]
		if len(values) > 1:
			row = {k: unidecode(v).strip() for k, v in zip(headers, values)}
			row['Sample'] = sample_name
			parsed_row = _parse_html_row(row)
			converted_table.append(parsed_row)
	return converted_table


def _parse_coverage(sample_name: str, coverage: BeautifulSoup) -> TableType:
	coverage_table = list()
	rows = coverage.find_all('tr')
	if len(rows) == 0:
		return coverage_table
	column_names = [i.text for i in rows[1].find_all('th')]

	for index, tag in enumerate(rows[2:]):
		values = tag.find_all('td')

		if len(values) > 1:
			row = [('Sample', sample_name)] + [(k, v.get_text()) for k, v in zip(column_names, values)]
			row = OrderedDict(row)

			row['start'] = to_integer(row['start'])
			row['end'] = to_integer(row['end'])
			row['size'] = to_integer(row['size'])
			coverage_table.append(row)

	return coverage_table


def _parse_junctions(sample_name: str, junctions: BeautifulSoup) -> TableType:
	rows = junctions.find_all('tr')
	column_names_a = ['0', '1'] + [unidecode(i.get_text()) for i in rows.pop(0).find_all('th')][1:]

	column_names_a[4] = '{} ({})'.format(column_names_a[4], 'single')
	column_names_b = [i for i in column_names_a if i not in ['reads (cov)', 'score', 'skew', 'freq', '0']]
	junction_table = list()
	for a_row, b_row in zip(rows[::2], rows[1::2]):
		a_values = [unidecode(i.get_text()) for i in a_row.find_all('td')]
		b_values = [unidecode(i.get_text()) for i in b_row.find_all('td')]

		a_row = {unidecode(k): v for k, v in zip(column_names_a, a_values)}
		b_row = {unidecode(k): v for k, v in zip(column_names_b, b_values)}
		a_row['Sample'] = sample_name
		b_row['Sample'] = sample_name
		junction_table.append(a_row)
		junction_table.append(b_row)
	return junction_table


def get_index_filename(path: Path) -> Path:
	path = Path(path)
	if path.suffix == '.html':
		return path
	index_file = path / "output" / "index.html"

	if not index_file.exists():
		candidates = list(path.glob("**/index.html"))
		if not candidates:
			message = f"Cannot find the index file for folder {path}"
			raise FileNotFoundError(message)
		elif len(candidates) != 1:
			message = f"Found multiple index files for folder {path}"
			raise FileNotFoundError(message)
		index_file = candidates[0]
	return index_file


def parse_index_file(sample_name: str, filename: Union[str, Path], set_index: bool = True, default_seq = 'chrom1') -> Tuple[DFType, DFType, DFType]:
	"""
		Extracts information on each of the tables from the index file.
	Parameters
	----------
	sample_name:
	filename
	set_index:bool; default True
		Whether to set the index of the dataframe.
	default_seq: str; default 'chrom1'
		Name of the sequence if the `seq id` column is not included in the output.

	Returns
	-------
	pandas.DataFrame
	- Index -> ()
	- Values-> ()
	"""
	filename = get_index_filename(filename)
	file_contents = _load_index_file(filename)

	snp_table_headers, snp_soup, coverage_soup, junction_soup = _extract_index_tables(file_contents)

	snp_table = _parse_snp_table(sample_name, snp_table_headers, snp_soup)

	try:
		coverage_table = _parse_coverage(sample_name, coverage_soup)
	# If the snp_table does not include a sequence id column, get it from the coverage table.
	except (ValueError, TypeError): coverage_table = []

	try:
		junction_table = _parse_junctions(sample_name, junction_soup)
	except (ValueError, TypeError, IndexError):
		junction_table = []

	snp_df = convert_to_dataframe(snp_table)
	coverage_df = convert_to_dataframe(coverage_table)
	junction_df = convert_to_dataframe(junction_table)

	if 'seq id' not in snp_df.columns:
		try:
			_sequence_id_from_coverage_table = coverage_df.iloc[0]['seq id']
		except (IndexError, KeyError):
			_sequence_id_from_coverage_table = default_seq
		snp_df['seq id'] = _sequence_id_from_coverage_table
	# Remove columns that shouldn't be there
	for col in snp_df.columns:
		if col not in VariantTableColumnMap:
			snp_df.pop(col)
	print(list(snp_df.columns))
	snp_df.columns = [VariantTableColumnMap[i] for i in snp_df.columns]
	if set_index:
		# Make sure the position column is a number. Breseq sometimes uses :1 if there is more than one mutation at a position.
		snp_df['position'] = [float(str(i).replace(':', '.').replace(',', '')) for i in snp_df['position']]
		snp_df.set_index(keys = ['seq id', 'position'], inplace = True)
	return snp_df, coverage_df, junction_df


