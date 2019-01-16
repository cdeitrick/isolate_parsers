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
from typing import NamedTuple, Optional, Tuple

import pandas

from breseqparser.file_parsers import parse_gd, parse_index, parse_vcf

DF = pandas.DataFrame
GDColumns = parse_gd.GDColumns
VCFColumns = parse_vcf.VCFColumns
IndexColumns = parse_index.VariantTableColumns


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
	alt: str = VCFColumns.alternate
	ref: str = VCFColumns.reference
	# alt:str = 'alt'
	# ref:str = 'ref'
	# quality: str = 'quality'
	# depth: str = 'readDepth'
	# variant_type: str = 'variantType'
	alternate_amino: str = GDColumns.alternate_amino
	reference_amino: str = GDColumns.reference_amino
	alternate_codon: str = GDColumns.alternate_codon
	reference_codon: str = GDColumns.reference_codon
	locus_tag: str = 'locusTag'
	mutation_category: str = 'mutationCategory'


IsolateTableColumns = _IsolateTableColumns()


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

def _filter_bp(raw_df: pandas.DataFrame) -> pandas.DataFrame:
	""" Filters out variants that occur within 1000bp of each other."""
	forward: pandas.Series = raw_df[IsolateTableColumns.position].diff().abs()
	reverse: pandas.Series = raw_df[IsolateTableColumns.position][::-1].diff()[::-1].abs()

	# noinspection PyTypeChecker
	fdf: pandas.DataFrame = raw_df[(forward > 1000) & (reverse > 1000) | (forward.isna() | reverse.isna())]

	return fdf

def parse_breseq_isolate(breseq_folder: Path, isolate_id: str, isolate_name: str = None, use_filter:bool = False) -> Tuple[DF, DF, DF]:
	"""
		Combines all available information for a single breseq (single sample).
	Parameters
	----------
	breseq_folder: Path
		Location of the breseq output folder. Should contain an index.html file, gd files, and a vcf file.
	isolate_id
	isolate_name

	Returns
	-------

	"""
	print("parsing ", isolate_id, "\t",breseq_folder)
	if not isolate_name:
		isolate_name = isolate_id

	_gd_subset_columns = [
		GDColumns.alternate_amino, GDColumns.reference_amino, GDColumns.alternate_codon,
		GDColumns.reference_codon, GDColumns.locus_tag, GDColumns.mutation_category,
		# GDColumns.alternate_base, GDColumns.reference_base
		# , GDColumns.position, GDColumns.sequence_id These columns are in the index.
	]

	# Retrieve the filenames
	index_file = parse_index.get_index_filename(breseq_folder)
	# The VCF file provides the quality and read depth.
	vcf_file = parse_vcf.get_vcf_filename(breseq_folder)
	gd_file = parse_gd.get_gd_filename(breseq_folder)

	# Generate the tables
	# They should be indexed by sequence id and position.
	_set_table_index = True
	index_df, coverage_df, junction_df = parse_index.parse_index_file(isolate_name, index_file, set_index = _set_table_index)
	vcf_df: pandas.DataFrame = parse_vcf.parse_vcf_file(vcf_file, set_index = _set_table_index, no_filter = True)
	gd_df = parse_gd.parse_gd_file(gd_file, set_index = _set_table_index)

	gd_subset = gd_df[_gd_subset_columns]

	variant_df: pandas.DataFrame = index_df.merge(gd_subset, how = 'left', left_index = True, right_index = True)
	variant_df = variant_df.merge(vcf_df, how = 'left', left_index = True, right_index = True)

	assert VCFColumns.alternate in variant_df.columns
	variant_df.reset_index(inplace = True)

	variant_df[IsolateTableColumns.sample_id] = isolate_id
	variant_df[IsolateTableColumns.sample_name] = isolate_name
	variant_df = variant_df[list(IsolateTableColumns)]

	if use_filter:
		variant_df = _filter_bp(variant_df)

	variant_df.set_index(keys = [IsolateTableColumns.sequence_id, IsolateTableColumns.position])

	return variant_df, coverage_df, junction_df


if __name__ == "__main__":
	path = Path(__file__).parent.parent / 'tests' / 'data' / 'Clonal_Output' / 'breseq_output'
	v, c, j = parse_breseq_isolate(path, 'testIsolate')
	print(v.to_string())
