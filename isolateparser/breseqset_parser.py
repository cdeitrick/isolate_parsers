""" Combines all subfolders generated from a set of breseq runs into a single spreadsheet."""

from pathlib import Path
from typing import Container, List, Mapping

import pandas
from loguru import logger

from isolateparser.breseqparser import BreseqOutputParser, get_sample_name



class BreseqIsolateSetParser:
	""" Parses a series of breseq calls for multiple isolates."""

	def __init__(self, whitelist: Container[str] = None, blacklist: Container[str] = None, sample_map: Mapping[str, str] = None,
			sample_regex: str = None, use_filter: bool = False):
		self.use_filter = use_filter
		self.whitelist = whitelist if whitelist else []
		self.blacklist = blacklist if whitelist else []
		self.sample_map = sample_map if sample_map else {}
		self.sample_regex = sample_regex  # To extract sample names from folders.

		self.variant_tables = list()
		self.coverage_tables = list()
		self.junction_table = list()
		self.summaries = list()

	@staticmethod
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
				if not f.exists():
					f = subfolder / "breseq"
				breseq_folders.append(f)

		return breseq_folders

	def run(self, parent_folder: Path):
		""" Expects a folder of breseq runs for a set ofisolates.
			Parameters
			----------
			parent_folder:Path
				The folder of breseq output folders.
		"""
		breseq_folders = self._get_breseq_folder_paths(parent_folder)

		for folder in breseq_folders:
			self.update_tables(folder)

		snp_dataframe_full = pandas.concat(self.variant_tables, sort = True)
		coverage_dataframe_full = pandas.concat(self.coverage_tables, sort = True)
		junction_dataframe_full = pandas.concat(self.junction_table, sort = True)

		summary = pandas.DataFrame(self.summaries)
		return snp_dataframe_full, coverage_dataframe_full, junction_dataframe_full, summary

	def update_tables(self, folder: Path):
		isolate_id = get_sample_name(folder)
		isolate_name = self.sample_map.get(isolate_id, isolate_id)
		in_whitelist = not self.whitelist or isolate_name in self.whitelist or isolate_id in self.whitelist
		in_blacklist = bool(self.blacklist) and (isolate_name in self.blacklist or isolate_id in self.blacklist)

		if in_blacklist or not in_whitelist: return None

		try:
			breseq_output = BreseqOutputParser(self.use_filter)
			snp_df, coverage_df, junction_df = breseq_output.run(folder, isolate_id, isolate_name)
			summary = breseq_output.get_summary(folder, isolate_id, isolate_name)
			self.variant_tables.append(snp_df)
			self.coverage_tables.append(coverage_df)
			self.junction_table.append(junction_df)
			self.summaries.append(summary)
		except FileNotFoundError as _missing_file_error:
			logger.warning(f"{_missing_file_error}")
			return None


