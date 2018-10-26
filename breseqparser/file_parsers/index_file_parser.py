from collections import OrderedDict
from pathlib import Path
from typing import List, Tuple, Any, Dict, Union
from pprint import pprint
from bs4 import BeautifulSoup
from unidecode import unidecode
import pandas
TableType = List[Dict[str,Any]]
DFType = pandas.DataFrame
def to_number(string: str) -> int:
	""" Converts a string to a number"""
	try:
		string = int(string.replace(',', ''))
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
	soup = BeautifulSoup(filename.read_text(), 'lxml')
	return soup


def _extract_index_file_tables(soup: BeautifulSoup) -> Tuple[List[str], BeautifulSoup, BeautifulSoup]:
	"""
		Extracts the headers for the snp table, the junction table, and the coverage table.
	Parameters
	----------
	soup: BeautifulSoup
		Equivilent to BeautifulSoup(index.html)

	Returns
	-------
		snp_header, coverage_table, junction_table

	"""
	alph_soup = str(soup)
	begin_snp_header_string = r'<th>evidence</th>'
	end_snp_header_string = '<!-- Item Lines -->'

	begin_snp_header = alph_soup.find(begin_snp_header_string)
	end_snp_header = alph_soup.find(end_snp_header_string)

	snp_header_full = alph_soup[begin_snp_header:end_snp_header]
	# snp_header = snp_header_full[end_snp_header:]

	snp_header_soup = BeautifulSoup(snp_header_full, 'lxml')
	snp_header_soup = [i.text for i in snp_header_soup.find_all('th')]

	begin_umc = alph_soup.find(
		'<tr><th align="left" class="missing_coverage_header_row" colspan="11">Unassigned missing coverage evidence</th></tr>')
	end_umc = alph_soup.find(
		'<th align="left" class="new_junction_header_row" colspan="12">Unassigned new junction evidence</th>')
	coverage_string = alph_soup[begin_umc:end_umc]
	junction_string = alph_soup[end_umc:]
	coverage_soup = BeautifulSoup(coverage_string, 'lxml')
	junction_soup = BeautifulSoup(junction_string, 'lxml')
	return snp_header_soup, coverage_soup, junction_soup


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

	normal_table = soup.find_all(attrs = {'class': 'normal_table_row'})
	poly_table = soup.find_all(attrs = {'class': 'polymorphism_table_row'})

	snp_table = normal_table + poly_table
	snp_header_soup, coverage_soup, junction_soup = _extract_index_file_tables(soup)

	return snp_header_soup, snp_table, coverage_soup, junction_soup


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
			row = {k: unidecode(v) for k, v in zip(headers, values)}
			try:
				row['Sample'] = sample_name
			except KeyError:
				row['Sample'] = "noname"
			try:
				row['position'] = to_number(row['position'])
			except KeyError:
				row['position'] = None
			try:
				row['freq %'] = float(row['freq'][:-1])
				row.pop('freq')
			except KeyError:
				pass
			if 'javascript' in row['description'] or 'Javascript' in row['description']:
				row['description'] = 'large deletion'
			converted_table.append(row)
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

			row['start'] = to_number(row['start'])
			row['end'] = to_number(row['end'])
			row['size'] = to_number(row['size'])
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

def get_index_filename(path:Path)->Path:
	path = Path(path)
	if path.name == 'index.html':
		return path
	index_file = path / "data" / "output" / "index.html"

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

def parse_index_file(sample_name: str, filename: Union[str,Path])->Tuple[DFType,DFType,DFType]:
	"""
		Extracts information on each of the tables from the index file.
	Parameters
	----------
	sample_name
	filename

	Returns
	-------

	"""
	filename = get_index_filename(filename)
	file_contents = _load_index_file(filename)

	snp_table_headers, snp_soup, coverage_soup, junction_soup = _extract_index_tables(file_contents)

	snp_table = _parse_snp_table(sample_name, snp_table_headers, snp_soup)
	coverage_table = _parse_coverage(sample_name, coverage_soup)
	junction_table = _parse_junctions(sample_name, junction_soup)

	snp_df = convert_to_dataframe(snp_table)
	coverage_df = convert_to_dataframe(coverage_table)
	junction_df = convert_to_dataframe(junction_table)

	return snp_df, coverage_df, junction_df

if __name__ == "__main__":
	path ="../data/index.html"

	_snp, _cov, _jun = parse_index_file("AU0074", path)
	print(_snp.to_string())

