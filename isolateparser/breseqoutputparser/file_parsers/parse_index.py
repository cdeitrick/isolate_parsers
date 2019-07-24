import math
from collections import OrderedDict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union

import pandas
from bs4 import BeautifulSoup
from unidecode import unidecode

TableType = List[Dict[str, Any]]
DFType = pandas.DataFrame

VariantTableColumnMap = {
	'Sample':      'sample',
	'annotation':  'annotation',
	'description': 'description',
	'freq %':      'frequency',
	'gene':        'gene',
	'mutation':    'mutation',
	'position':    'position',
	'seq id':      'seq id'
}


# General Utilities

def to_integer(string: str) -> int:
	""" Converts a string to a number"""
	try:
		string = int(float(string.replace(',', '')))
	except ValueError:
		pass
	return string


def _read_index_file(filename: Path) -> BeautifulSoup:
	"""Loads the index file."""
	filename = Path(filename)
	# Use read_bytes to avoid encoding errors.
	soup = BeautifulSoup(filename.read_bytes(), 'lxml')
	return soup


###################################################################################################################################################
########################################################## General Table Parsing Methods ##########################################################
###################################################################################################################################################
def _fix_html_fieldnames(columns: Iterable[str]) -> List[str]:
	""" The `seq id` field is encoded weirdly due to html, so need to fix it. """
	return [(i if i != 'seq\xa0id' else 'seq id') for i in columns]


###################################################################################################################################################
################################################################ SNP Table parsing ################################################################
###################################################################################################################################################
def _extract_variant_table_headers(text: str) -> List[str]:
	begin_snp_header_string = r'<th>evidence</th>'
	end_snp_header_string = '<!-- Item Lines -->'

	begin_snp_header = text.find(begin_snp_header_string)
	end_snp_header = text.find(end_snp_header_string)

	snp_header_full = text[begin_snp_header:end_snp_header]

	snp_header_soup = BeautifulSoup(snp_header_full, 'lxml')
	snp_header_soup = [i.text for i in snp_header_soup.find_all('th')]
	return snp_header_soup


def _parse_variant_table_html_row(row: Dict[str, str]) -> Dict[str, str]:
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


def _parse_snp_table(headers: List[str], rows: List[BeautifulSoup]) -> TableType:
	"""
		Parses the SNP table.
	Parameters
	----------
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
			parsed_row = _parse_variant_table_html_row(row)
			converted_table.append(parsed_row)
	return converted_table


def _extract_snp_table(soup: BeautifulSoup) -> List[BeautifulSoup]:
	normal_table = soup.find_all(attrs = {'class': 'normal_table_row'})
	poly_table = soup.find_all(attrs = {'class': 'polymorphism_table_row'})
	snp_table = normal_table + poly_table

	return snp_table


def add_missing_columns(snp_df: pandas.DataFrame, coverage_df: pandas.DataFrame, default_seq: str, sample_name) -> pandas.DataFrame:
	# Sometimes `seq id` is not in the snp table, so it needs to be added manually for compatibility with other scripts.
	if 'seq id' not in snp_df.columns:
		try:
			_sequence_id_from_coverage_table = coverage_df.iloc[0]['seq id']
		except (IndexError, KeyError):
			_sequence_id_from_coverage_table = default_seq
		snp_df['seq id'] = _sequence_id_from_coverage_table

	if 'freq' not in snp_df: snp_df['freq'] = math.nan  # Should unpack into a sequence automatically

	snp_df['sampleName'] = sample_name

	return snp_df


###################################################################################################################################################
####################################################### Coverage and Junction Table Parsing #######################################################
###################################################################################################################################################

def _extract_coverage_and_junction_tables(alph_soup: str) -> Tuple[BeautifulSoup, BeautifulSoup]:
	"""
		Extracts the headers for the snp table, the junction table, and the coverage table.
	Parameters
	----------
	alph_soup: str

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


def _parse_coverage_table_row(tag: BeautifulSoup, column_names: List[str]) -> Optional[Dict[str, str]]:
	""" Converts a single row in the coverage table into a dictionary."""
	values = tag.find_all('td')

	if len(values) > 1:
		row = [(k, v.get_text()) for k, v in zip(column_names, values)]
		row = OrderedDict(row)

		row['start'] = to_integer(row['start'])
		row['end'] = to_integer(row['end'])
		row['size'] = to_integer(row['size'])
		return row


def _parse_coverage_soup(coverage: BeautifulSoup) -> TableType:
	rows = coverage.find_all('tr')
	if len(rows) == 0:
		return []
	column_names = [i.text for i in rows[1].find_all('th')]
	coverage_table = list()
	for index, tag in enumerate(rows[2:]):
		row = _parse_coverage_table_row(tag, column_names)
		if row:
			coverage_table.append(row)

	return coverage_table


def _parse_junction_table_row_pair(first: BeautifulSoup, second: BeautifulSoup, names_first, names_second) -> Tuple[Dict[str, str], Dict[str, str]]:
	a_values = [unidecode(i.get_text()) for i in first.find_all('td')]
	b_values = [unidecode(i.get_text()) for i in second.find_all('td')]

	a_row = {unidecode(k): v for k, v in zip(names_first, a_values)}
	b_row = {unidecode(k): v for k, v in zip(names_second, b_values)}

	return a_row, b_row


def _parse_junction_soup(junctions: BeautifulSoup) -> TableType:
	rows = junctions.find_all('tr')
	if not len(rows): return []

	# Extract the column names from the junction table.
	column_names_a = ['0', '1'] + [unidecode(i.get_text()) for i in rows.pop(0).find_all('th')][1:]
	column_names_a[4] = '{} ({})'.format(column_names_a[4], 'single')
	column_names_b = [i for i in column_names_a if i not in ['reads (cov)', 'score', 'skew', 'freq', '0']]

	junction_table = list()
	for a_row, b_row in zip(rows[::2], rows[1::2]):
		parsed_row_a, parsed_row_b = _parse_junction_table_row_pair(a_row, b_row, column_names_a, column_names_b)
		junction_table += [parsed_row_a, parsed_row_b]

	return junction_table


###################################################################################################################################################
################################################################ Main Parser ######################################################################
###################################################################################################################################################
def parse_index_file(sample_name: str, filename: Union[str, Path], set_index: bool = True, default_seq = 'chrom1') -> Tuple[DFType, DFType, DFType]:
	"""
		Extracts information on each of the tables from the index file.
	Parameters
	----------
	sample_name:
	filename: Path
		Path to the index.html file. Assume this links to the file itself and delegate file finding to another process.
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
	file_soup = _read_index_file(filename)
	file_contents = str(file_soup)

	snp_table_headers = _extract_variant_table_headers(file_contents)
	snp_soup = _extract_snp_table(file_soup)
	coverage_soup, junction_soup = _extract_coverage_and_junction_tables(file_contents)

	snp_table = _parse_snp_table(snp_table_headers, snp_soup)

	try:
		coverage_table = _parse_coverage_soup(coverage_soup)
	# If the snp_table does not include a sequence id column, get it from the coverage table.
	except (ValueError, TypeError): coverage_table = []

	try:
		junction_table = _parse_junction_soup(junction_soup)
	except (ValueError, TypeError, IndexError):
		junction_table = []

	# convert each table to a DataFrame
	snp_df = pandas.DataFrame(snp_table)
	coverage_df = pandas.DataFrame(coverage_table)
	junction_df = pandas.DataFrame(junction_table)

	# The `seq id` fieldname is encoded weirdly in the html file, so need ot fix it.
	snp_df.columns = _fix_html_fieldnames(snp_df.columns)
	coverage_df.columns = _fix_html_fieldnames(coverage_df.columns)
	junction_df.columns = _fix_html_fieldnames(junction_df.columns)

	# Remap the column names to something a little more readable.
	snp_df.columns = [VariantTableColumnMap.get(i, i) for i in snp_df.columns]

	# Add missing columns for compatibility.
	snp_df = add_missing_columns(snp_df, coverage_df, default_seq, sample_name)
	coverage_df['sampleName'] = sample_name
	junction_df['sampleName'] = sample_name

	if set_index:
		# Make sure the position column is a number. Breseq sometimes uses :1 if there is more than one mutation at a position.
		snp_df['position'] = [float(str(i).replace(':', '.').replace(',', '')) for i in snp_df['position']]
		snp_df.set_index(keys = ['seq id', 'position'], inplace = True)

	return snp_df, coverage_df, junction_df


class IndexFileParser:
	def __init__(self, sample_name: str):
		self.sample_name = sample_name
