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
from loguru import logger

from .parsers import GDParser, IndexParser, parse_summary_file, parse_vcf_file

DF = pandas.DataFrame


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
	alt: str = 'alt'
	ref: str = 'ref'
	alternate_amino: str = 'aminoAlt'
	reference_amino: str = 'aminoRef'
	alternate_codon: str = 'codonAlt'
	reference_codon: str = 'codonRef'
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


class BreseqOutputParser:
	def __init__(self, use_filter: bool = False):
		self._set_table_index = True
		self.use_filter = use_filter

		self.index_columns = []
		# We only want a subset of the columns available from the gd  and vcf files.
		self.gd_columns = ['aminoAlt', 'aminoRef', 'codonAlt', 'codonRef', 'mutationCategory']
		self.vcf_columns = ['alt', 'ref']

		self.file_parser_index = IndexParser()
		self.file_parser_gd = GDParser()

	def run(self, sample_id: str, indexpath: Path, gdpath: Optional[Path] = None, vcfpath: Optional[Path] = None,
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

		index_df, coverage_df, junction_df = self.file_parser_index.run(sample_name, indexpath, set_index = self._set_table_index)

		if gdpath:
			gd_df = self.file_parser_gd.run(gdpath, set_index = self._set_table_index)
		else:
			# TODO: Need to add the missing columns
			# These would be 'locusTag' and 'mutationCategory'.
			gd_df = None

		if vcfpath:
			vcf_df = parse_vcf_file(vcfpath, set_index = self._set_table_index, no_filter = True)
		else:
			# TODO: Need to add the missing columns
			# These would be the 'alt' and 'ref' columns
			vcf_df = None

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
			# We don't care about most of the columns
			reduced_gd = gd[self.gd_columns]

			variant_df: pandas.DataFrame = index.merge(reduced_gd, how = 'left', left_index = True, right_index = True)
		else:
			variant_df = index
		logger.warning(f"Saving the df for debugging")
		variant_df.to_csv("varianttesttable.tsv", sep = '\t')
		if vcf is not None:
			variant_df = variant_df.merge(vcf, how = 'left', left_index = True, right_index = True)

		return variant_df

	@staticmethod
	def get_summary(folder: Path, sample_id: str, sample_name: Optional[str] = None) -> Dict[str, Any]:
		return parse_summary_file(folder, sample_id, sample_name)

	def get_alt_from_gd(self):
		pass
