from pathlib import Path
from typing import List, Dict
import pandas

from breseqparsers.isolate_parser import parse_breseq_isolate


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


def save_isolate_table(tables:Dict[str,pandas.DataFrame], filename: Path) -> Path:
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
		writer = pandas.ExcelWriter(filename)
		include_index = False
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

def parse_breseqset(folder: Path):
	""" Expects a folder of breseq runs for a set ofisolates."""
	breseq_folders = _get_breseq_folder_paths(folder)
	snp_dfs = list()
	coverage_dfs = list()
	junction_dfs = list()
	for index, breseq_folder in enumerate(breseq_folders):
		try:
			snp_df, coverage_df, junction_df = parse_breseq_isolate(breseq_folder)
		except FileNotFoundError:
			continue
		snp_dfs.append(snp_df)
		coverage_dfs.append(coverage_df)
		junction_dfs.append(junction_df)


	snp_dataframe_full = pandas.concat(snp_dfs)
	coverage_dataframe_full = pandas.concat(coverage_dfs)
	junction_dataframe_full = pandas.concat(junction_dfs)

	tables = {
		'snp': snp_dataframe_full,
		'coverage': coverage_dataframe_full,
		'junction': junction_dataframe_full
	}
	save_isolate_table(tables, folder / "breseq_table.xlsx")

	# snp_dataframe_full.to_csv(folder / "breseq_table.tsv", sep = "\t")

if __name__ == "__main__":
	_folder = Path("/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/lipuma/pipeline_output/")

	parse_breseqset(_folder)
