import argparse
from pathlib import Path
from typing import *
import pandas
from loguru import logger


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


	if arguments:
		_program_options = parser.parse_args(arguments)
	else:
		_program_options = parser.parse_args()
	return _program_options
