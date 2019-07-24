import argparse
from pathlib import Path
from typing import Dict, List, Union

from loguru import logger


def load_program_options(arguments: List[str] = None) -> argparse.Namespace:
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"-i", "--input",
		help = "The breseq folder to parse.",
		dest = "folder",
		type = Path
	)
	parser.add_argument(
		"--no-fasta",
		help = "Whether to generate an aligned fasta file of all snps in the breseq VCF file.",
		action = 'store_false',
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
	if arguments:
		_program_options = parser.parse_args(arguments)
	else:
		_program_options = parser.parse_args()
	return _program_options


class IsolateSetWorkflow:
	__version__ = "0.1.0"

	def __init__(self, arguments: argparse.Namespace):
		self.whitelist = self._parse_commandline_list(arguments.whitelist)
		self.blacklist = self._parse_commandline_list(arguments.blacklist)
		self.fasta_categories = self._parse_commandline_list(arguments.snp_categories)
		self.use_filter = arguments.use_filter
		self.generate_fasta = arguments.generate_fasta

		if arguments.sample_map: self.sample_map = self._parse_sample_map(arguments.sample_map)
		else: self.sample_map = {}

		self.breseq_callset_parser = BreseqCallSetParser(
			whitelist = self.whitelist,
			blacklist = self.blacklist,
			sample_map = self.sample_map,
			use_filter = self.use_filter
		)

	def run(self, parent_folder: Path, reference_label: str):
		prefix = parent_folder.name
		output_filename_table = parent_folder / f"{prefix}.xlsx"
		output_filename_fasta = parent_folder / f"{prefix}"

		variant_df, coverage_df, junction_df, summary_df = self.breseq_callset_parser.run(parent_folder)

		logger.info("Generating comparison table...")
		snp_comparison_df = generate_snp_comparison_table(variant_df, 'base', self.use_filter, reference_label)
		amino_comparison_df = generate_snp_comparison_table(variant_df, 'amino', self.use_filter, reference_label)
		codon_comparison_df = generate_snp_comparison_table(variant_df, 'codon', self.use_filter, reference_label)

		tables = {
			'variant comparison': snp_comparison_df,
			'amino comparison':   amino_comparison_df,
			'codon comparison':   codon_comparison_df,
			'variant':            variant_df.reset_index(),
			'coverage':           coverage_df.reset_index(),
			'junction':           junction_df.reset_index(),
			'summary':            summary_df
		}
		logger.info("Saving isolate table as ", output_filename_table)
		save_isolate_table(tables, output_filename_table)

		if self.generate_fasta:
			logger.info("Generating fasta...")
			fasta_filename_snp = output_filename_fasta.with_suffix(".snp.fasta")
			generate_fasta_file(variant_df, fasta_filename_snp, by = 'base', reference_label = program_options.reference_label)

	@staticmethod
	def _parse_commandline_list(io: Union[None, str, List[str]]) -> List[str]:
		"""
			Attempts to convert a comma-separated list of options given from the command line.
		Parameters
		----------
		io: str
			Either a comma-separated list of ids or a file path of a text file with each id occupying a single line.
		Returns
		-------
		List[str]
		"""
		if isinstance(io, list): return io
		filename = Path(io)

		# Make sure io is not ''
		if io and filename.exists():
			contents = filename.read_text().split('\n')
		else:
			# Assume it is a comma-separated list.
			# An empty string will result in an empty list.
			contents = io.split(',')
		contents = [i for i in contents if i]
		return contents

	@staticmethod
	def _parse_sample_map(path: str) -> Dict[str, str]:
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
			print(message)
			contents = {}

		if 'sampleId' in contents:
			# The column name isn't necessary, but won't break the parser.
			# pop the key since it isn't needed.
			contents.pop('sampleId')
		return contents


if __name__ == "__main__":
	from isolateparser.breseq_callset_parser import BreseqCallSetParser
	from isolateparser.generate import generate_snp_comparison_table, save_isolate_table, generate_fasta_file

	debug_options = [
		"--input", "/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/lipuma/pipeline/SC1360/",
		"--sample-map", "/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/lipuma/isolate_sample_map.old.tsv",
		"--reference", "A0-01"
	]

	program_options = load_program_options()
	isolateset_workflow = IsolateSetWorkflow(program_options)
	isolateset_workflow.run(program_options.folder, program_options.reference_label)
