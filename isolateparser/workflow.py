from pathlib import Path
from loguru import logger
from typing import *
from isolateparser.resultparser import BreseqFolderParser, get_sample_name
from isolateparser.resultparser.parsers import locations
from isolateparser.generate import generate_fasta_file, generate_snp_comparison_table, save_isolate_table
import pandas
import argparse

class IsolateSetWorkflow:
	"""
		Parses a folder containing numerious breseq call folders. The directory structure should
		look like this:
		parent_folder
			sampleId
				breseq
					data
					output
			sampleId
				breseq
					data
					output
	"""
	__version__ = "0.1.0"

	def __init__(self, whitelist: Union[str, List[str]] = "", blacklist: Union[str, List[str]] = "", sample_map: Path = None,
			sample_regex: str = None,
			use_filter: bool = False, snp_categories: Container[str] = None, generate_fasta: bool = True):
		self.whitelist = self._parse_commandline_list(whitelist)
		if blacklist is None:
			self.blacklist = []
		else:
			self.blacklist = self._parse_commandline_list(blacklist)
		self.fasta_categories = self._parse_commandline_list(snp_categories)

		self.use_filter = use_filter
		self.generate_fasta = generate_fasta
		self.sample_regex = sample_regex

		# This decides which categories of mutations are used to generate the aligned fasta file.
		self.snp_categories = snp_categories if snp_categories else ["snp_synonymous", "snp_nonsynonymous"]

		if sample_map: self.sample_map = self._parse_sample_map(sample_map)
		else: self.sample_map = {}

		self.variant_tables = list()
		self.coverage_tables = list()
		self.junction_table = list()
		self.gd_extra_tables = list()
		self.summaries = list()

	def run(self, parent_folder: Path, reference_label: str, output_prefix: Optional[Path] = None) -> Path:
		if output_prefix:
			output_prefix = output_prefix.absolute()
			prefix = output_prefix.name
			output_filename_table = output_prefix.parent / f"{prefix}.xlsx"
			output_filename_fasta = output_prefix.parent / f"{prefix}"
		else:
			prefix = parent_folder.absolute().name
			output_filename_table = parent_folder / f"{prefix}.xlsx"
			output_filename_fasta = parent_folder / f"{prefix}"

		variant_df, coverage_df, junction_df, summary_df = self.concatenate_callset_tables(parent_folder)

		logger.info("Generating comparison table...")
		snp_comparison_df = generate_snp_comparison_table(variant_df, 'base', reference_label)

		tables = {
			'variant comparison': snp_comparison_df,
			'variant':            variant_df.reset_index(),
			'coverage':           coverage_df.reset_index(),
			'junction':           junction_df.reset_index(),
			'summary':            summary_df
		}
		# Omit `None` from the list.
		self.gd_extra_tables = [i for i in self.gd_extra_tables if i is not None]
		if self.gd_extra_tables:
			gd_df = pandas.concat(self.gd_extra_tables, sort = False)
			tables['gd'] = gd_df

		logger.info(f"Saving isolate table as {output_filename_table}")
		save_isolate_table(tables, output_filename_table)

		if self.generate_fasta:
			logger.info("Generating fasta...")
			is_population = variant_df['frequency'].nunique() > 1
			if is_population:
				message = "Generating fastas is not supported for populations and may fail."
				logger.warning(message)

			fasta_filename_snp = output_filename_fasta.with_suffix(".snp.fasta")
			fasta_filename_codon = output_filename_fasta.with_suffix(".codon.fasta")
			generate_fasta_file(variant_df, fasta_filename_snp, by = 'base', reference_label = reference_label)
			generate_fasta_file(variant_df, fasta_filename_codon, by = 'codon', reference_label = reference_label)

		return output_filename_table

	def concatenate_callset_tables(self, parent_folder: Path):
		""" Expects a folder of breseq runs for a set of isolates.
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
		if self.summaries:
			summary = pandas.DataFrame(self.summaries)
		else:
			summary = None
		return snp_dataframe_full, coverage_dataframe_full, junction_dataframe_full, summary

	@staticmethod
	def _get_breseq_folder_paths(base_folder: Path) -> List[Path]:
		""" Attempts to find all folders corresponding to a breseq run.
		Situation 1
		.parent
		|---- sample 1
		|----|---- breseq
		|----|----|---- output

		Situation 2
		.parent
		|---- sample 1
		|----|---- output

		Situation 3
		.parent
		|---- output

		"""
		breseq_folders = list()
		for subfolder in base_folder.iterdir():
			if subfolder.is_file(): continue
			result = locations.is_folder_breseq(subfolder)
			# The index.html file is the bare minimum.
			if result:
				breseq_folders.append(subfolder)
			else:
				message = f"Cannot find an index.html file in {subfolder}. Skipping..."
				logger.warning(message)
		return breseq_folders

	@staticmethod
	def _parse_commandline_list(data: Union[None, str, List[str]]) -> List[str]:
		"""
			Attempts to convert a comma-separated list of options given from the command line.
		Parameters
		----------
		data: str
			Either a comma-separated list of ids or a file path of a text file with each id occupying a single line.
		Returns
		-------
		List[str]
		"""
		if data is None: data = []
		if isinstance(data, list): return data
		filename = Path(data)

		# Make sure io is not ''
		if data and filename.exists():
			contents = filename.read_text().split('\n')
		else:
			# Assume it is a comma-separated list.
			# An empty string will result in an empty list.
			contents = data.split(',')
		contents = [i for i in contents if i]
		return contents

	@staticmethod
	def _parse_sample_map(path: Union[None, str, Path]) -> Dict[str, str]:
		"""
			Converts a file mapping sample ids to sample names into a usable dictionary.
		Parameters
		----------
		path: str
			path to the file.
		Returns
		-------
		Dict[str,str]
		Maps sample ids to sample names.
		"""
		if not path: return {}
		try:
			filename = Path(path)

			contents = dict()
			lines = filename.read_text().split('\n')
			for line in lines:
				# Check for extra empty lines.
				if not line: continue
				key, value = line.split('\t')
				contents[key] = value
		except FileNotFoundError:
			contents = {}
		except ValueError:
			message = "The sample map file is not formatted correctly. Make sure all lines contain exactly two values: the sample id and sample name."
			logger.warning(message)
			contents = {}

		if 'sampleId' in contents:
			# The column name isn't necessary, but won't break the parser.
			# pop the key since it isn't needed.
			contents.pop('sampleId')
		return contents

	def update_tables(self, folder: Path):
		isolate_id = get_sample_name(folder)
		isolate_name = self.sample_map.get(isolate_id, isolate_id)
		logger.info(f"Parsing '{isolate_name}' ('{isolate_id}')")
		in_whitelist = not self.whitelist or isolate_name in self.whitelist or isolate_id in self.whitelist
		in_blacklist = bool(self.blacklist) and (isolate_name in self.blacklist or isolate_id in self.blacklist)

		if in_blacklist or not in_whitelist: return None

		breseq_output = BreseqFolderParser(self.use_filter)

		filenames_breseq = locations.get_file_locations(folder)
		verbose = False
		if verbose:
			logger.info(f"Files in folder {folder}")
			logger.info(f"\tIndex file: {filenames_breseq['index']}")
			logger.info(f"\tVcf file: {filenames_breseq['vcf']}")
			logger.info(f"\tGd file: {filenames_breseq['gd']}")
			logger.info(f"\tSummary: {filenames_breseq['summary']}")

		snp_df, coverage_df, junction_df, gd_extra_df = breseq_output.run(
			indexpath = filenames_breseq['index'],
			gdpath = filenames_breseq['gd'],
			vcfpath = filenames_breseq['vcf'],
			sample_id = isolate_id,
			sample_name = isolate_name
		)

		self.variant_tables.append(snp_df)
		self.coverage_tables.append(coverage_df)
		self.junction_table.append(junction_df)
		if gd_extra_df is not None:
			gd_extra_df['sampleId'] = isolate_id
			gd_extra_df['sampleName'] = isolate_name
			self.gd_extra_tables.append(gd_extra_df)

		if filenames_breseq.get('summary'):
			summary = breseq_output.get_summary(filenames_breseq['summary'], isolate_id, isolate_name)
			self.summaries.append(summary)


class IsolateParser:
	# TODO: This class isn't being used for anything.
	""" A parser that is dedicated to parsing a single folder."""

	def __init__(self, sample_map: Dict[str, str] = None):
		self.use_filter = False
		self.sample_map = sample_map if sample_map else {}

	def run(self, folder, output_filename: Path = None):
		isolate_id = get_sample_name(folder)
		isolate_name = self.sample_map.get(isolate_id, isolate_id)
		breseq_output = BreseqFolderParser(self.use_filter)

		filenames_breseq = locations.get_file_locations(folder)

		logger.info(f"Files in folder {folder}")
		logger.info(f"\tIndex file: {filenames_breseq['index']}")
		logger.info(f"\tVcf file: {filenames_breseq['vcf']}")
		logger.info(f"\tGd file: {filenames_breseq['gd']}")
		logger.info(f"\tSummary: {filenames_breseq['summary']}")

		snp_df, coverage_df, junction_df = breseq_output.run(
			indexpath = filenames_breseq['index'],
			gdpath = filenames_breseq['gd'],
			vcfpath = filenames_breseq['vcf'],
			sample_id = isolate_id,
			sample_name = isolate_name
		)
		logger.debug(f"Length of coverage: {len(coverage_df)}")
		logger.debug(f"Length pf Junction: {len(junction_df)}")

		summary = breseq_output.get_summary(folder, isolate_id, isolate_name)
		summary_df = pandas.Series(summary).to_frame().reset_index().transpose()

		tables = {
			'variant':  snp_df.reset_index(),
			'coverage': coverage_df.reset_index(),
			'junction': junction_df.reset_index(),
			'summary':  summary_df
		}
		for key, value in tables.items():
			logger.debug(f"{key}: {type(value)}")
		if output_filename is None:
			output_filename = "breseq.xlsx"
		logger.info(f"Saving isolate table as {output_filename}")
		save_isolate_table(tables, output_filename)

def create_parser(arguments: List[str] = None) -> argparse.Namespace:
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"-i", "--input",
		help = "The breseq folder to parse.",
		dest = "folder",
		type = Path
	)
	parser.add_argument(
		"-o", "--output",
		help = "Where to save the output files. Should just be the prefix, the file extensions will be added automatically.",
		dest = "output",
		type = Path,
		default = None
	)
	parser.add_argument(
		"--fasta",
		help = "Whether to generate an aligned fasta file of all snps in the breseq VCF file.",
		action = 'store_true',
		dest = 'generate_fasta'
	)

	parser.add_argument(
		"-w", "--whitelist",
		help = "Samples not in the whitelist are ignored. Either a comma-separated list of sample ids for a file with each sample id occupying a single line.",
		dest = "whitelist",
		action = 'store',
		default = ""
	)

	parser.add_argument(
		"-b", "--blacklist",
		help = "Samples to ignore. See `--whitelist` for possible input formats.",
		action = 'store',
		dest = 'blacklist',
		default = ""
	)
	parser.add_argument(
		"-m", "--sample-map",
		help = """A file mapping sample ids to sample names. Use if the subfolders in the breseqset folder are named differently from the sample names."""
			   """ The file should have two columns: `sampleId` and `sampleName`, separated by a tab character.""",
		action = 'store',
		dest = 'sample_map',
		default = ""
	)
	parser.add_argument(
		"--filter-1000bp",
		help = "Whether to filter out variants that occur within 1000bp of each other. Usually indicates a mapping error.",
		action = "store_true",
		dest = "use_filter"
	)
	parser.add_argument(
		"--reference",
		help = "The sample that was used as the reference, if available.",
		action = "store",
		default = None,
		dest = "reference_label"
	)
	parser.add_argument("--snp-categories", help = "Categories to use when concatenating SNPs into a fasta file.", dest = "snp_categories",
		default = "")

	parser.add_argument("--regex", help = "Used to extract sample names from the given filename. Currently Disabled", type = str)

	parser.add_argument("--single", help = "Indicates that there is only one sample. Used for debugging.", action = "store_true")

	if arguments:
		_program_options = parser.parse_args(arguments)
	else:
		_program_options = parser.parse_args()

	return _program_options

def main():
	program_options = create_parser()

	if program_options.single:
		isolate_workflow = IsolateParser(program_options.sample_map)
		isolate_workflow.run(
			folder = program_options.folder.absolute(),
			output_filename = program_options.output
		)
	else:
		isolateset_workflow = IsolateSetWorkflow(
			whitelist = program_options.whitelist,
			blacklist = program_options.blacklist,
			sample_map = program_options.sample_map,
			sample_regex = program_options.regex,
			use_filter = program_options.use_filter,
			snp_categories = program_options.snp_categories,
			generate_fasta = program_options.generate_fasta
		)
		isolateset_workflow.run(program_options.folder, program_options.reference_label, program_options.output)
if __name__ == "__main__":
	main()
