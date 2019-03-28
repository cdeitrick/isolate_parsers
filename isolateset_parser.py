import argparse
from pathlib import Path
from typing import Dict, List, Optional, Union

from dataclasses import dataclass


@dataclass
class ProgramOptions:
	folder: str
	generate_fasta: bool = True
	whitelist: List[str] = ""
	blacklist: List[str] = ""
	sample_map: str = ""
	use_filter: bool = True
	reference_label: Optional[str] = None


def _get_program_options(arguments:List[str] = None) -> Union[ProgramOptions, argparse.Namespace]:
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"-i", "--input",
		help = "The breseq folder to parse.",
		dest = "folder"
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
	if arguments:
		_program_options = parser.parse_args(arguments)
	else:
		_program_options = parser.parse_args()
	return _program_options


def _parse_commandline_list(io: Optional[str]) -> List[str]:
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
	from breseqset_parser import parse_breseqset
	from file_generators import generate_snp_comparison_table, save_isolate_table, generate_fasta_file, generate_basic_statistics

	program_options = _get_program_options()
	whitelist = _parse_commandline_list(program_options.whitelist)
	blacklist = _parse_commandline_list(program_options.blacklist)
	sample_map = _parse_sample_map(program_options.sample_map)

	breseq_run_folder = Path(program_options.folder)
	breseq_table_filename = breseq_run_folder / "breseq_table.xlsx"
	fasta_filename_base = breseq_run_folder / "breseq"
	print("Parsing breseqset...")
	variant_df, coverage_df, junction_df = parse_breseqset(breseq_run_folder, blacklist, whitelist, sample_map, program_options.use_filter)
	assert 'ref' in variant_df
	print("Generating comparison table...")
	snp_comparison_df = generate_snp_comparison_table(variant_df, by = 'base', filter_table = program_options.use_filter)
	amino_comparison_df = generate_snp_comparison_table(variant_df, by = 'amino', filter_table = program_options.use_filter)
	codon_comparison_df = generate_snp_comparison_table(variant_df, by = 'codon', filter_table = program_options.use_filter)
	assert 'ref' in variant_df

	tables = {
		'variant comparison': snp_comparison_df,
		'amino comparison':   amino_comparison_df,
		'codon comparison':   codon_comparison_df,
		'variant':            variant_df.reset_index(),
		'coverage':           coverage_df.reset_index(),
		'junction':           junction_df.reset_index()
	}
	print("Saving isolate table...")
	save_isolate_table(tables, breseq_run_folder / "breseq_table.xlsx")

	if program_options.generate_fasta:
		print("Generating fasta...")
		fasta_filename_snp = fasta_filename_base.with_suffix(f".snp.fasta")
		fasta_filename_codon = fasta_filename_base.with_suffix(f".codon.fasta")
		fasta_filename_amino = fasta_filename_base.with_suffix(f".amino.fasta")
		generate_fasta_file(variant_df, fasta_filename_snp, by = 'base', reference_label = program_options.reference_label)
		generate_fasta_file(variant_df, fasta_filename_codon, by = 'codon', reference_label = program_options.reference_label)
		generate_fasta_file(variant_df, fasta_filename_amino, by = 'amino', reference_label = program_options.reference_label)
