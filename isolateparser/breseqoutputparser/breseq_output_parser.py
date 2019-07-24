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
from typing import Any, Dict, NamedTuple, Optional, Tuple

import pandas

from isolateparser.breseqoutputparser.file_parsers import locations, parse_gd, parse_index, parse_summary, parse_vcf

DF = pandas.DataFrame
GDColumns = parse_gd.GDColumns
VCFColumns = parse_vcf.VCFColumns


class _IsolateTableColumns(NamedTuple):
	# Defines the column names for the isolate table. Used `Enum` because it's iterable and only one instance will exist.
	sample_id: str = 'sampleId'
	sample_name: str = 'sampleName'
	sequence_id: str = 'seq id'
	position: str = 'position'
	annotation: str = 'annotation'
	description: str = 'description'
	evidence: str = 'evidence'
	freq: str = 'freq'
	gene: str = 'gene'
	mutation: str = 'mutation'
	alt: str = VCFColumns.alternate
	ref: str = VCFColumns.reference
	alternate_amino: str = GDColumns.alternate_amino
	reference_amino: str = GDColumns.reference_amino
	alternate_codon: str = GDColumns.alternate_codon
	reference_codon: str = GDColumns.reference_codon
	locus_tag: str = 'locusTag'
	mutation_category: str = 'mutationCategory'


IsolateTableColumns = _IsolateTableColumns()


def get_sample_name(folder: Path) -> Optional[str]:
	""" Attempt to extract the sample name from a folder."""
	return [i for i in folder.parts if 'breseq' not in i][-1]


def string_is_date(value: str) -> bool:
	length = len(value)
	result = length > 8 and len([i for i in value if i.isdigit()]) > (len(value) / 2)
	return result


def _filter_bp(raw_df: pandas.DataFrame) -> pandas.DataFrame:
	# TODO: This should only be done for sequences on the same chrom.
	""" Filters out variants that occur within 1000bp of each other."""
	forward: pandas.Series = raw_df[IsolateTableColumns.position].diff().abs()
	reverse: pandas.Series = raw_df[IsolateTableColumns.position][::-1].diff()[::-1].abs()

	# noinspection PyTypeChecker
	fdf: pandas.DataFrame = raw_df[(forward > 1000) & (reverse > 1000) | (forward.isna() | reverse.isna())]

	return fdf


def get_file_locations(folder: Path) -> Tuple[Path, Optional[Path], Optional[Path], Optional[Path]]:
	index_file = locations.get_index_filename(folder)
	gd_file = locations.get_gd_filename(folder)
	# The VCF file provides the quality and read depth.
	vcf_file = locations.get_vcf_filename(folder)

	summary_file = parse_summary.get_summary_filename(folder)

	return index_file, vcf_file, gd_file, summary_file


class BreseqOutputParser:
	def __init__(self, use_filter: bool = False):
		self._set_table_index = True
		self.use_filter = use_filter

		self.index_columns = []
		# We only want a subset of the columns available from the gd file.
		self.gd_columns = [
			GDColumns.alternate_amino, GDColumns.reference_amino, GDColumns.alternate_codon,
			GDColumns.reference_codon, GDColumns.locus_tag, GDColumns.mutation_category,
		]
		self.vcf_columns = []

	def run(self, sample_id: str, indexpath: Path, vcfpath: Optional[Path] = None, gdpath: Optional[Path] = None,
			sample_name: Optional[str] = None) -> Tuple[
		pandas.DataFrame, pandas.DataFrame, pandas.DataFrame]:
		"""
			Runs the workflow.
			Parameters
			----------
				sample_id:str
				indexpath: Path
				vcfpath: Optional[Path]
				gdpath: Optional[Path]
				sample_name: Optional[str]
		"""
		if not sample_name:
			sample_name = sample_id

		index_df, coverage_df, junction_df = parse_index.parse_index_file(sample_name, indexpath, set_index = self._set_table_index)
		if vcfpath:
			vcf_df = parse_vcf.parse_vcf_file(vcfpath, set_index = self._set_table_index, no_filter = True)
		else:
			vcf_df = None
		if gdpath:
			gd_df = parse_gd.parse_gd_file(gdpath, set_index = self._set_table_index)
		else:
			gd_df = None
		# Merge the tables together.
		variant_df = self.merge_tables(index_df, gd_df, vcf_df)

		# Add the `sampleId` and `sampleName`
		variant_df['sampleId'] = sample_id
		variant_df['sampleName'] = sample_name
		variant_df = variant_df[[i for i in IsolateTableColumns if i in variant_df.columns]]

		return variant_df, coverage_df, junction_df

	def merge_tables(self, index: pandas.DataFrame, gd: Optional[pandas.DataFrame], vcf: Optional[pandas.DataFrame]) -> pandas.DataFrame:
		"""
			Merges the three tables that contain mutational information.
		Parameters
		----------
		index: pandas.DataFrane
			The dataframe representing the index.html file.
		gd: Optional[pandas.DataFrame]
			The converted form of the gd file, if available. If the data cannot be found, this will not be merged.
		vcf: Optional[pandsa.DataFrame]
			Same as for the gd file.

		Returns
		-------
		pandas.DataFrame
			A dataframe that contains data from all three given tables.
		"""
		if gd is not None:
			variant_df: pandas.DataFrame = index.merge(gd[self.gd_columns], how = 'left', left_index = True, right_index = True)
		else:
			variant_df = index

		if vcf is not None:
			variant_df = variant_df.merge(vcf, how = 'left', left_index = True, right_index = True)

		return variant_df

	@staticmethod
	def get_summary(folder: Path, sample_id: str, sample_name: Optional[str] = None) -> Dict[str, Any]:
		return parse_summary.parse_summary_file(folder, sample_id, sample_name)
