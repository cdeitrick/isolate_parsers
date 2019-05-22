from pathlib import Path
from typing import Any, Dict, List, NamedTuple

import pandas
import vcf


class _VCFColumns(NamedTuple):
	sequence_id: str = 'seq id'
	position: str = 'position'
	alternate: str = 'alt'
	reference: str = 'ref'
	quality: str = 'quality'
	depth: str = 'readDepth'
	variant_type: str = 'variantType'


VCFColumns = _VCFColumns()


def _filter_df(raw_df: pandas.DataFrame) -> pandas.DataFrame:
	""" Filters out variants that occur within 1000bp of each other."""
	forward: pandas.Series = raw_df[VCFColumns.position].diff().abs()
	reverse: pandas.Series = raw_df[VCFColumns.position][::-1].diff()[::-1].abs()

	# noinspection PyTypeChecker
	fdf: pandas.DataFrame = raw_df[(forward > 1000) & (reverse > 1000) | (forward.isna() | reverse.isna())]

	return fdf


def get_vcf_filename(path: Path) -> Path:
	""" Attempts to locate the vcf file in the outputfolder generated by breseq.
		Parameters
		----------
		path: Path
			The folder generated by breseq.
		Returns
		-------
		Path
			The vcf file generated by breseq.
	"""
	if path.is_dir():
		result = path / "data" / "output.vcf"
		if not result.exists():
			candidates = list(path.glob("**/output.vcf"))
			if len(candidates) != 1:
				message = f"Invalid vcf file path: {path}"
				raise FileNotFoundError(message)
			result = candidates[0]
	elif path.suffix == '.vcf':
		result = path
	else:
		message = f"Invalid VCF path: {path}"
		raise ValueError(message)
	return result


def _convert_record_to_dictionary(record: Any) -> Dict[str, str]:
	alt = "".join(str(i) for i in record.ALT)
	ref = "".join(str(i) for i in record.REF)
	qual = record.QUAL
	depth = record.INFO.get('DP')
	record_data = {
		VCFColumns.sequence_id:  record.CHROM,
		VCFColumns.position:     record.POS,
		VCFColumns.alternate:    alt,
		VCFColumns.reference:    ref,
		VCFColumns.quality:      qual,
		VCFColumns.depth:        depth,
		VCFColumns.variant_type: record.var_type
	}

	return record_data


def _convert_vcf_to_table(vcf_filename: Path) -> List[Dict[str, Any]]:
	"""Converts all records in a vcf file into a list of dictionaries."""
	table: List[Dict[str, str]] = list()
	seen_positions = set()
	with vcf_filename.open('r') as file1:
		vcf_reader = vcf.Reader(file1)
		for record in vcf_reader:
			data = _convert_record_to_dictionary(record)
			#VCF files sometimes record separate mutations as occuring at the same position.
			# The gd file will instead increment the second mutations position by 1. Do this to maintain compatibility.
			if (data['seq iq'], data['position']) in seen_positions:
				data['position'] += 1
			seen_positions.add((data['seq id'], data['position']))
			table.append(data)
	return table


def parse_vcf_file(path: Path, set_index: bool = True, no_filter: bool = False) -> pandas.DataFrame:
	"""
		Converts the VCF file generated by breseq into a pandas Dataframe.
	Parameters
	----------
	path: Path
		Either a folder containing a single breseq run or a path to the vcf file itself.
	set_index:bool; default True
		Whether to set the index of the dataframe.
	no_filter: bool default False
		Whether to use the filters or not. By default, variants occuring within 1000bp of each other are excluded.
	Returns
	-------
	pandas.DataFrame
		- Index -> (VCFColumns.sample_name, VCFColumns.sequence_id, VCFColumns.position)
		- Columns-> VCFColumns
	"""
	filename = get_vcf_filename(path)

	table = _convert_vcf_to_table(filename)

	# Columns are defined in VCFColumns
	vcf_df: pandas.DataFrame = pandas.DataFrame(table)

	# Filter out variants that occur within 1000b bo of each other.
	if no_filter:
		filtered_df = vcf_df
	else:
		filtered_df = _filter_df(vcf_df)

	# Order the columns correctly
	filtered_df = filtered_df[list(VCFColumns)]

	if set_index:
		filtered_df.set_index(keys = [VCFColumns.sequence_id, VCFColumns.position], inplace = True)

	return filtered_df
