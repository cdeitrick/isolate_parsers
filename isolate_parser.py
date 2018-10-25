"""
	Parses the results from a single breseq run agains a single isolate.
	Expects the following structure:
	./breseq_output/
	├── 01_sequence_conversion
	├── 02_reference_alignment
	├── 03_candidate_junctions
	├── 04_candidate_junction_alignment
	├── 05_alignment_correction
	├── 06_bam
	├── 07_error_calibration
	├── 08_mutation_identification
	├── data
	│   ├── AU0074.forward.trimmed.paired.unmatched.fastq
	│   ├── AU0074.forward.trimmed.unpaired.unmatched.fastq
	│   ├── AU0074.reverse.trimmed.paired.unmatched.fastq
	│   ├── output.gd
	│   ├── output.vcf
	│   ├── reference.bam
	│   ├── reference.bam.bai
	│   ├── reference.fasta
	│   ├── reference.fasta.fai
	│   ├── reference.gff3
	│   └── summary.json
	├── output
	│   ├── calibration
	│   ├── evidence
	│   ├── index.html
	│   ├── log.txt
	│   ├── marginal.html
	│   ├── output.done
	│   ├── output.gd
	│   └── summary.html
	└── sample_output

"""

from pathlib import Path
from index_file_parser import parse_index_file
from typing import Optional
import re
from pprint import pprint
def _get_index_file_path(path:Path)->Path:
	if path.name == 'index.html':
		return path
	index_file = path / "data" / "output" / "index.html"

	if not index_file.exists():
		candidates = list(path.glob("**/index.html"))
		if len(candidates) != 1:
			message = f"Cannot find the index file for folder {path}"
			pprint(candidates)
			raise FileNotFoundError(message)
		index_file = candidates[0]
	return index_file


def _get_sample_name(folder:Path)->Optional[str]:
	""" Attempt to extract the sample name from a folder."""
	pattern = "[A-Z]+[\d]+"
	for part in folder.parts[::-1]:
		match = re.search(pattern, part)
		if match:
			name = part
			break
	else:
		name = None

	return name
def parse_breseq_isolate(breseq_folder:Path):
	""" parses all available data from a single breseq run."""
	isolate_name = _get_sample_name(breseq_folder)
	index_file = _get_index_file_path(breseq_folder)

	snp_df, coverage_df, junction_df = parse_index_file(isolate_name, index_file)

	snp_df['sampleName'] = coverage_df['sampleName'] = junction_df['sampleName'] = isolate_name
	return snp_df, coverage_df, junction_df



