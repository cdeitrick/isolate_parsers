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

import re
from pathlib import Path
from typing import Optional

from breseqparser.file_parsers import vcf_file_parser, index_file_parser


def _get_sample_name(folder: Path) -> Optional[str]:
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


def parse_breseq_isolate(breseq_folder: Path):
	""" parses all available data from a single breseq run."""
	isolate_name = _get_sample_name(breseq_folder)
	index_file = index_file_parser.get_index_filename(breseq_folder)
	vcf_file = vcf_file_parser.get_vcf_filename(breseq_folder)

	snp_df, coverage_df, junction_df = index_file_parser.parse_index_file(isolate_name, index_file)
	vcf_df = vcf_file_parser.parse_vcf(vcf_file, isolate_name)

	snp_df['sampleName'] = coverage_df['sampleName'] = junction_df['sampleName'] = isolate_name
	snp_df = snp_df.set_index(keys = ['sampleName', 'seq id', 'position'])

	variant_df = snp_df.merge(vcf_df, how = 'left', left_index = True, right_index = True)
	return variant_df, coverage_df, junction_df
