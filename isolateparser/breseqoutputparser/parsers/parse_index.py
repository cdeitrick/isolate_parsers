import math
import re
from collections import OrderedDict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union

import pandas
from bs4 import BeautifulSoup
from loguru import logger
from unidecode import unidecode

TableType = List[Dict[str, Any]]
DFType = pandas.DataFrame


###################################################################################################################################################
########################################################## General Table Parsing Methods ##########################################################
###################################################################################################################################################
def _clean_fieldnames(columns: Iterable[str]) -> List[str]:
	""" The `seq id` field is encoded weirdly due to html, so need to fix it. """
	return [(i if i != 'seq\xa0id' else 'seq id') for i in columns]


def to_integer(string: str) -> int:
	""" Converts a string to a number"""
	try:
		string = int(float(string.replace(',', '')))
	except ValueError:
		pass
	return string


###################################################################################################################################################
################################################################ Variant Table parsing ############################################################
###################################################################################################################################################


class VariantTableParser:
	column_map = {
		'Sample':      'sample',
		'annotation':  'annotation',
		'description': 'description',
		'freq %':      'frequency',
		'gene':        'gene',
		'position':    'position',
		'seq id':      'seq id'
	}

	def __init__(self):
		self.header_string_start = r'<th>evidence</th>'
		self.header_string_end = '<!-- Item Lines -->'

	@staticmethod
	def _combine_dictionary_values(left: Dict[str, str], right: Dict[str, str]) -> Dict[str, str]:
		available_keys = sorted(set(left.keys()) | set(right.keys()))
		combined = dict()
		for key in available_keys:
			left_value = left.get(key)
			right_value = right.get(key)

			if left_value and not right_value:
				value = left_value
			elif right_value and not left_value:
				value = right_value
			elif left_value and right_value:
				if left_value == right_value:
					value = left_value
				else:
					value = f"{left_value}|{right_value}"
			else: value = None
			combined[key] = value
		return combined

	def _extract_fields_from_translated_cds_annotation(self, annotation: str) -> Dict[str, str]:
		# [locus_tag=Bcen2424_5147] [db_xref=InterPro:IPR006157] [protein=Dihydroneopterin aldolase][protein_id=ABK11880.1] [location=complement(2194414..2194854)][gbkey=CDS]/[locus_tag=Bcen2424_5146] [db_xref=InterPro:IPR008258] [protein=Lytic transglycosylase,catalytic] [protein_id=ABK11879.1] [location=2193627..2194313][gbkey=CDS]
		if ']/[' in annotation:
			# This annotation is refering to multiple locus tags
			left, right = self._split_annotation(annotation)
			left_result = self._extract_fields_from_translated_cds_annotation(left)
			right_result = self._extract_fields_from_translated_cds_annotation(right)

			return self._combine_dictionary_values(left_result, right_result)

		tokenize_pattern = "\[([a-z_]+)[=](.+?)\]"
		matches = re.findall(tokenize_pattern, annotation)

		return dict(matches)

	@staticmethod
	def _extract_snp_table(soup: BeautifulSoup) -> List[BeautifulSoup]:
		normal_table = soup.find_all(attrs = {'class': 'normal_table_row'})
		poly_table = soup.find_all(attrs = {'class': 'polymorphism_table_row'})
		snp_table = normal_table + poly_table

		return snp_table

	def _extract_table_headers(self, soup: Union[str, BeautifulSoup]) -> List[str]:
		text = str(soup)  # do not use prettify()
		begin_snp_header = text.find(self.header_string_start)
		end_snp_header = text.find(self.header_string_end)

		snp_header_full = text[begin_snp_header:end_snp_header]

		snp_header_soup = BeautifulSoup(snp_header_full, 'lxml')
		snp_header_soup = [i.text for i in snp_header_soup.find_all('th')]
		return snp_header_soup

	def _parse_snp_table(self, headers: List[str], rows: List[BeautifulSoup]) -> TableType:
		"""
			Parses the SNP table.
		Parameters
		----------
		headers: List[str]
			Column names for the snp table.
		rows: List
			The rows of the snp table.
		"""
		converted_table = list()
		for tag in rows:
			values = [v.text for v in tag.find_all('td')]
			if len(values) > 1:
				row = {k: unidecode(v).strip() for k, v in zip(headers, values)}
				parsed_row = self._parse_variant_table_html_row(row)
				converted_table.append(parsed_row)
		return converted_table

	def _parse_variant_table_html_row(self, row: Dict[str, Union[str, int, float]]) -> Dict[str, str]:
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
			row['position'] = math.nan
		try:
			row['freq %'] = float(row['freq'][:-1])
		# row.pop('freq')
		except KeyError:
			pass
		# `description` is the column name in the index file.

		if 'javascript' in row['description'] or 'Javascript' in row['description']:
			row['description'] = 'large_deletion'
		if row['description'] == 'large_deletion' or row['description'] == 'large deletion':
			row['alt'] = 'deletion'

		#logger.info(row)

		# Now add annotations from the `description` field if it was annotated from a translated cds file.
		# [locus_tag=Bcen2424_5147] [db_xref=InterPro:IPR006157] [protein=Dihydroneopterin aldolase][protein_id=ABK11880.1] [location=complement(2194414..2194854)][gbkey=CDS]/[locus_tag=Bcen2424_5146] [db_xref=InterPro:IPR008258] [protein=Lytic transglycosylase,catalytic] [protein_id=ABK11879.1] [location=2193627..2194313][gbkey=CDS]

		description = row['description']
		if 'locus_tag' in description:
			description = row['description']
			fields = self._extract_fields_from_translated_cds_annotation(description)

			if 'locus_tag' in fields:
				row['locusTag'] = fields['locus_tag']

			row = {**row, **fields}

		return row

	@staticmethod
	def _split_annotation(annotation: str) -> Tuple[str, str]:
		""" Splits a cds annotation that refers to a mutaiton that spans multiple locus tags."""
		# Assume that ']/[' is the separator

		left, right = annotation.split(']/[')
		left = left + ']'
		right = '[' + right
		return left, right

	def run(self, sample_name: str, soup: BeautifulSoup) -> pandas.DataFrame:
		snp_table_headers = self._extract_table_headers(str(soup))
		snp_table = self._extract_snp_table(soup)

		parsed_table = self._parse_snp_table(snp_table_headers, snp_table)

		snp_df = pandas.DataFrame(parsed_table)

		snp_df.columns = _clean_fieldnames(snp_df.columns)
		# make sure the `position` column is type int
		snp_df['position'] = snp_df['position'].apply(lambda s: s if ':' not in s else s.split(':')[0])
		snp_df['position'] = snp_df['position'].astype(int)

		# Remap the column names to something a little more readable.
		snp_df.columns = [self.column_map.get(i, i) for i in snp_df.columns]

		snp_df['sampleName'] = sample_name
		return snp_df


###################################################################################################################################################
####################################################### Coverage and Junction Table Parsing #######################################################
###################################################################################################################################################

class ExtraTableParser:
	def __init__(self):
		self.umc_table_start_string = '<tr><th align="left" class="missing_coverage_header_row" colspan="11">Unassigned missing coverage evidence</th></tr>'
		self.umc_table_end_string = '<th align="left" class="new_junction_header_row" colspan="12">Unassigned new junction evidence</th>'

	def run(self, sample_name: str, soup: BeautifulSoup) -> Tuple[pandas.DataFrame, pandas.DataFrame]:
		contents = soup.prettify()

		coverage_soup, junction_soup = self._extract_coverage_and_junction_tables(contents)
		coverage_table = self._parse_coverage_soup(coverage_soup)
		junction_table = self._parse_junction_soup(junction_soup)
		coverage_df = pandas.DataFrame(coverage_table)
		junction_df = pandas.DataFrame(junction_table)

		coverage_df.columns = _clean_fieldnames(coverage_df.columns)
		junction_df.columns = _clean_fieldnames(junction_df.columns)

		coverage_df['sampleName'] = sample_name
		junction_df['sampleName'] = sample_name

		return coverage_df, junction_df

	def _extract_coverage_and_junction_tables(self, alph_soup: str) -> Tuple[BeautifulSoup, BeautifulSoup]:
		"""
			Extracts the headers for the snp table, the junction table, and the coverage table.
		Parameters
		----------
		alph_soup: str

		Returns
		-------
			coverage_table, junction_table
		"""

		umc_table_index_start = alph_soup.find(
			self.umc_table_start_string)
		umc_table_index_end = alph_soup.find(
			self.umc_table_end_string)

		coverage_string = alph_soup[umc_table_index_start:umc_table_index_end]
		junction_string = alph_soup[umc_table_index_end:]
		coverage_soup = BeautifulSoup(coverage_string, 'lxml')
		junction_soup = BeautifulSoup(junction_string, 'lxml')
		return coverage_soup, junction_soup

	def _parse_coverage_soup(self, coverage: BeautifulSoup) -> TableType:
		rows = coverage.find_all('tr')
		if len(rows) == 0:
			return []
		column_names = [i.text for i in rows[1].find_all('th')]
		coverage_table = list()
		for index, tag in enumerate(rows[2:]):
			row = self._parse_coverage_table_row(tag, column_names)
			if row:
				coverage_table.append(row)

		return coverage_table

	@staticmethod
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

	def _parse_junction_soup(self, junctions: BeautifulSoup) -> TableType:
		rows = junctions.find_all('tr')
		if not len(rows): return []

		# Extract the column names from the junction table.
		column_names_a = ['0', '1'] + [unidecode(i.get_text()) for i in rows.pop(0).find_all('th')][1:]
		column_names_a[4] = '{} ({})'.format(column_names_a[4], 'single')
		column_names_b = [i for i in column_names_a if i not in ['reads (cov)', 'score', 'skew', 'freq', '0']]

		junction_table = list()
		for a_row, b_row in zip(rows[::2], rows[1::2]):
			parsed_row_a, parsed_row_b = self._parse_junction_table_row_pair(a_row, b_row, column_names_a, column_names_b)
			junction_table += [parsed_row_a, parsed_row_b]

		return junction_table

	@staticmethod
	def _parse_junction_table_row_pair(first: BeautifulSoup, second: BeautifulSoup, names_first, names_second) -> Tuple[
		Dict[str, str], Dict[str, str]]:
		a_values = [unidecode(i.get_text()) for i in first.find_all('td')]
		b_values = [unidecode(i.get_text()) for i in second.find_all('td')]

		a_row = {unidecode(k): v for k, v in zip(names_first, a_values)}
		b_row = {unidecode(k): v for k, v in zip(names_second, b_values)}

		return a_row, b_row


###################################################################################################################################################
################################################################ Main Parser ######################################################################
###################################################################################################################################################


class IndexParser:
	def __init__(self, default_seq: str = 'chrom1'):
		"""

		Parameters
		----------
		default_seq: str; default 'chrom1'
			Name of the sequence if the `seq id` column is not included in the output.
		"""
		self.variant_table_parser = VariantTableParser()
		self.cj_table_parser = ExtraTableParser()

		self.default_seq = default_seq  # Only used if the seq id is not already available.

	def run(self, sample_name: str, filename: Path, set_index: bool = True):
		"""
			Extracts information on each of the tables from the index file.
		Parameters
		----------
		sample_name:
		filename: Path
			Path to the index.html file. Assume this links to the file itself and delegate file finding to another process.
		set_index:bool; default True
			Whether to set the index of the dataframe.


		Returns
		-------
		pandas.DataFrame
		- Index -> ()
		- Values-> ()
		"""
		file_soup = self._read_index_file(filename)
		variant_table = self.variant_table_parser.run(sample_name, file_soup)
		coverage_table, junction_table = self.cj_table_parser.run(sample_name, file_soup)

		variant_table = self.add_missing_columns(variant_table, coverage_table, self.default_seq)

		if set_index:
			variant_table = self.set_index(variant_table)
		return variant_table, coverage_table, junction_table

	@staticmethod
	def _read_index_file(filename: Path) -> BeautifulSoup:
		"""Loads the index file."""
		filename = Path(filename)
		# Use read_bytes to avoid encoding errors.
		soup = BeautifulSoup(filename.read_bytes(), 'lxml')
		return soup

	@staticmethod
	def add_missing_columns(snp_df: pandas.DataFrame, coverage_df: pandas.DataFrame, default_seq: str) -> pandas.DataFrame:
		# Sometimes `seq id` is not in the snp table, so it needs to be added manually for compatibility with other scripts.
		if 'seq id' not in snp_df.columns:
			try:
				_sequence_id_from_coverage_table = coverage_df.iloc[0]['seq id']
			except (IndexError, KeyError):
				_sequence_id_from_coverage_table = default_seq
			snp_df['seq id'] = _sequence_id_from_coverage_table

		if 'freq' not in snp_df: snp_df['freq'] = math.nan  # Should unpack into a sequence automatically

		return snp_df

	@staticmethod
	def set_index(table: pandas.DataFrame) -> pandas.DataFrame:
		# Make sure the position column is a number. Breseq sometimes uses :1 if there is more than one mutation at a position.
		table['position'] = [float(str(i).replace(':', '.').replace(',', '')) for i in table['position']]
		table.set_index(keys = ['seq id', 'position'], inplace = True)
		return table
