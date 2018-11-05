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

from breseqparser.file_parsers import parse_vcf, parse_index
from breseqparser.file_parsers import parse_gd


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


def parse_breseq_isolate(breseq_folder: Path, isolate_name:str = None):
	""" parses all available data from a single breseq run."""
	if not isolate_name:
		isolate_name = _get_sample_name(breseq_folder)
	index_file = parse_index.get_index_filename(breseq_folder)
	vcf_file = parse_vcf.get_vcf_filename(breseq_folder)
	gd_file = parse_gd.get_gd_filename(breseq_folder)
	snp_df, coverage_df, junction_df = parse_index.parse_index_file(isolate_name, index_file)
	vcf_df = parse_vcf.parse_vcf(vcf_file, isolate_name)
	gd_df = parse_gd.parse_gd(gd_file, isolate_name)
	gd_subset = gd_df[['aminoAlt', 'aminoRef', 'locusTag', 'mutationCategory', 'position', 'seq id', 'sampleName']]
	gd_subset = gd_subset.set_index(keys = ['sampleName', 'seq id', 'position'])
	snp_df['sampleName'] = coverage_df['sampleName'] = junction_df['sampleName'] = isolate_name
	snp_df = snp_df.set_index(keys = ['sampleName', 'seq id', 'position'])

	variant_df = snp_df.merge(vcf_df, how = 'left', left_index = True, right_index = True)
	variant_df = variant_df.merge(gd_subset, how = 'left', left_index = True, right_index = True)
	return variant_df, coverage_df, junction_df

if __name__ == "__main__":
	path = Path(__file__).parent.parent / 'data' / "breseq_run" / "AU0074"
	v, *_ = parse_breseq_isolate(path, 'AU0074')
	print(v.to_string())
