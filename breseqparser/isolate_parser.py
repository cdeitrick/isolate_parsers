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
from typing import NamedTuple, Optional

import pandas

from breseqparser.file_parsers import parse_gd, parse_index, parse_vcf


class _IsolateTableColumns(NamedTuple):
	# Defines the column names for the isolate table. Used `Enum` because it's iterable and only one instance will exist.
	sample_id: str = 'sampleId'
	sample_name: str = 'sampleName'
	sequence_id: str = 'seq id'
	position: str = 'position'
	annotation: str = 'annotation'
	description: str = 'description'
	evidence: str = 'evidence'
	gene: str = 'gene'
	mutation: str = 'mutation'
	alt: str = 'alt'
	quality: str = 'quality'
	depth: str = 'readDepth'


IsolateTableColumns = _IsolateTableColumns()
GDColumns = parse_gd.GDColumns
VCFColumns = parse_vcf.VCFColumns
IndexColumns = parse_index.VariantTableColumns


def get_sample_name(folder: Path) -> Optional[str]:
	""" Attempt to extract the sample name from a folder."""
	use_pattern = False
	if use_pattern:
		pattern = "[A-Z]+[\d]+|WSPlus|WC"
		for part in folder.parts[::-1]:
			match = re.search(pattern, part)
			if match:
				name = part
				break
		else:
			name = None
	else:
		return [i for i in folder.parts if 'breseq' not in i][-1]
	return name


def parse_breseq_isolate(breseq_folder: Path, isolate_id: str, isolate_name: str = None):
	""" parses all available data from a single breseq run."""

	if not isolate_name:
		isolate_name = isolate_id

	_gd_subset_columns = [
		GDColumns.alternate_amino, GDColumns.reference_amino, GDColumns.locus_tag,
		GDColumns.mutation_category  # , GDColumns.position, GDColumns.sequence_id These columns are in the index.
	]

	# Retrieve the filenames
	index_file = parse_index.get_index_filename(breseq_folder)
	vcf_file = parse_vcf.get_vcf_filename(breseq_folder)
	gd_file = parse_gd.get_gd_filename(breseq_folder)

	# Generate the tables
	# They should be indexed by sequence id and position.
	_set_table_index = True
	index_df, coverage_df, junction_df = parse_index.parse_index_file(isolate_name, index_file, set_index = _set_table_index)
	vcf_df: pandas.DataFrame = parse_vcf.parse_vcf_file(vcf_file, set_index = _set_table_index, no_filter = True)
	gd_df = parse_gd.parse_gd_file(gd_file, set_index = _set_table_index)

	gd_subset = gd_df[_gd_subset_columns]

	variant_df = index_df.merge(vcf_df, how = 'left', left_index = True, right_index = True)
	variant_df = variant_df.merge(gd_subset, how = 'left', left_index = True, right_index = True)

	variant_df[IsolateTableColumns.sample_id] = isolate_id
	variant_df[IsolateTableColumns.sample_name] = isolate_name

	return variant_df, coverage_df, junction_df


if __name__ == "__main__":
	path = Path(__file__).parent.parent / 'data' / "breseq_run" / "AU0074"
	v, *_ = parse_breseq_isolate(path, 'AU0074')
	print(v.to_string())
