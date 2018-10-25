from pathlib import Path
from typing import List, Optional

import pandas

from isolate_parser import parse_breseq_isolate


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


def _get_index_file_path(breseq_folder: Path) -> Optional[Path]:
	index_file = breseq_folder / "data" / "output" / "index.html"

	if not index_file.exists():
		candidates = list(breseq_folder.glob("**/index.html"))
		if len(candidates) != 1:
			message = "Cannot find the index file for folder {}".format(breseq_folder)
			print(message)
			return None
		index_file = candidates[0]
	return index_file


def parse_breseqset(folder: Path):
	""" Expects a folder of breseq runs for a set ofisolates."""
	breseq_folders = _get_breseq_folder_paths(folder)
	snp_dfs = list()
	coverage_dfs = list()
	junction_dfs = list()
	for index, breseq_folder in enumerate(breseq_folders):
		index_filename = _get_index_file_path(breseq_folder)
		if index_filename is None:
			continue
		snp_df, coverage_df, junction_df = parse_breseq_isolate(index_filename)
		snp_dfs.append(snp_df)
		coverage_dfs.append(coverage_df)
		junction_dfs.append(junction_df)

	snp_dataframe_full = pandas.concat(snp_dfs)
	coverage_dataframe_full = pandas.concat(coverage_dfs)
	junction_dataframe_full = pandas.concat(junction_dfs)

	snp_dataframe_full.to_excel(folder / "breseq_table.xlsx")


	# snp_dataframe_full.to_csv(folder / "breseq_table.tsv", sep = "\t")

if __name__ == "__main__":
	_folder = Path("/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/lipuma/pipeline_output/")

	parse_breseqset(_folder)
