""" Combines all subfolders generated from a set of breseq runs into a single spreadsheet."""

from pathlib import Path
from typing import Container, Dict, List

import pandas

from breseqparser.isolate_parser import get_sample_name, parse_breseq_isolate
from file_generators import save_isolate_table


def _get_breseq_folder_paths(base_folder: Path) -> List[Path]:
	""" Attempts to find all folders corresponding to a breseq run."""
	breseq_folders = list()
	for subfolder in base_folder.iterdir():
		if subfolder.is_file(): continue
		folder_contents = list(i.name for i in subfolder.iterdir())
		if 'output' in folder_contents and 'data' in folder_contents:
			# The provided folder is a folder of breseq runs.
			breseq_folders.append(subfolder)
		else:
			# Assume it is a folder of sample folders, each containing a `breseq_output` folder.
			f = subfolder / "breseq_output"
			breseq_folders.append(f)

	return breseq_folders


def parse_breseqset(folder: Path, blacklist: Container[str] = None, whitelist: Container[str] = None, sample_map: Dict[str, str] = None, use_filter:bool = False):
	""" Expects a folder of breseq runs for a set ofisolates.
		Parameters
		----------
		folder:Path
			The folder of breseq output folders.
		blacklist: List[str]
			A list of sample ids to ignore when generating the table.
		whitelist: List[str]
			Sample Ids not in `whitelist` will be excluded.
		sample_map: Dict[str,str]
	"""

	if not blacklist: blacklist = []
	if not whitelist: whitelist = []
	if not sample_map: sample_map = {}

	breseq_folders = _get_breseq_folder_paths(folder)
	snp_dfs = list()
	coverage_dfs = list()
	junction_dfs = list()
	for breseq_folder in breseq_folders:
		isolate_id = get_sample_name(breseq_folder)
		isolate_name = sample_map.get(isolate_id, isolate_id)

		try:
			snp_df, coverage_df, junction_df = parse_breseq_isolate(
				breseq_folder,
				isolate_id = isolate_id,
				isolate_name = isolate_name,
				use_filter = use_filter
			)
		except FileNotFoundError as _missing_file_error:
			print(f"Exception: {_missing_file_error}")
			continue

		in_whitelist = not whitelist or isolate_name in whitelist or isolate_id in whitelist
		in_blacklist = bool(blacklist) and (isolate_name in blacklist or isolate_id in blacklist)

		if in_blacklist: continue
		elif not in_whitelist: continue

		snp_dfs.append(snp_df)
		coverage_dfs.append(coverage_df)
		junction_dfs.append(junction_df)

	snp_dataframe_full = pandas.concat(snp_dfs, sort = True)

	chromosomes = snp_dataframe_full.groupby(by = ['seq id', 'position'])
	mutation_counts = dict()
	for (chrom, position), group in chromosomes:
		mutation_counts[chrom, position] = len(group)
	# snp_dataframe_full['presentInNSamples'] = [	mutation_counts[j, k] for j, k in snp_dataframe_full.index]

	coverage_dataframe_full = pandas.concat(coverage_dfs, sort = True)
	junction_dataframe_full = pandas.concat(junction_dfs, sort = True)

	return snp_dataframe_full, coverage_dataframe_full, junction_dataframe_full


if __name__ == "__main__":
	_folder = Path("/media/cld100/FA86364B863608A1/Users/cld100/Storage/TravisanoBreseq/")
	_variant_dataframe, _coverage_dataframe, _junction_dataframe = parse_breseqset(_folder)

	_tables = {
		'snp':      _variant_dataframe.reset_index(),
		'coverage': _coverage_dataframe.reset_index(),
		'junction': _junction_dataframe.reset_index()
	}

	save_isolate_table(_tables, _folder / "breseq_table.xlsx")
