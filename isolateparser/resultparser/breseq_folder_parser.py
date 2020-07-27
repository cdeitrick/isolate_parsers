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
from typing import *
import pandas
from loguru import logger

from .parsers import GDParser, IndexParser, parse_summary_file, parse_vcf_file, GDToTable
from ..tableformat import IsolateTableColumns

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

def catagorize_mutation(annotation:str, mutation:str)->str:
	""" Used when the mutation category is not available from the gd file."""
	if '-' in annotation:
		# is a SNP
		if 'intergenic' in annotation:
			result = "snp_intergenic"
		else:
			# Should be an amino acid annotation
			# Ex. G117G (GGC→GGG)
			aminoannotation = annotation.split(' ')[0]
			result = "snp_synonymous" if aminoannotation[0] == aminoannotation[-1] else "snp_nonsynonymous"
	elif annotation.startswith('(') or annotation.startswith('+'):
		# Should be an insertion.
		# (GGTGCCC)1
		# (GGTGCCC)1→2
		result = "small_indel"
	elif annotation.startswith('D'):
		# Δ425,094 bp
		# Try to check the size of the deletion
		numbers = [i for i in annotation if i.isdigit()]
		length = int("".join(numbers))
		if length < 100:
			result = 'deletion'
		else:
			result = "large_deletion"

	else:
		logger.debug(f"Could not categorize '{mutation}', '{annotation}'")
		result = "unknown"

	return result

def get_reference_from_mutation(mutation:str, annotation:str)->str:
	if '-' in mutation and len(mutation) == 3:
		ref = mutation.split('-')[0]
	elif mutation.startswith('+'):
		# An insertion
		ref = "."
	elif mutation.startswith('('):
		# Should be an insertion
		ref = "".join([i for i in mutation if i in "ACGT"])
		logger.warning(ref)
	elif annotation.startswith('D'):
		# It's a deletion
		ref = "."
	else:
		ref = '.'

	return ref

def get_alternate_from_mutation(mutation:str)->str:
	unknown_value = "N/A"
	if '-' in mutation and len(mutation) == 3:
		result = mutation.split('-')[-1]
	elif mutation.startswith('+'):
		result = mutation[1]
	elif mutation.startswith('D'):
		result = mutation
	else:
		result = unknown_value
	return result


class BreseqFolderParser:
	"""
		Parses the contents of a breseq run folder.
	"""
	def __init__(self, use_filter: bool = False):
		self._set_table_index = True
		self.use_filter = use_filter

		self.index_columns = []
		# We only want a subset of the columns available from the gd  and vcf files.
		self.gd_columns = ['aminoAlt', 'aminoRef', 'codonAlt', 'codonRef', 'mutationCategory']
		self.vcf_columns = ['alt', 'ref', 'readDepth']

		self.file_parser_index = IndexParser()
		self.file_parser_gd = GDParser()
		self.file_parser_gd_extra = GDToTable()

	def run(self, sample_id: str, indexpath: Path, gdpath: Optional[Path] = None, vcfpath: Optional[Path] = None,
			sample_name: Optional[str] = None) -> Tuple[
		pandas.DataFrame, pandas.DataFrame, pandas.DataFrame, pandas.DataFrame]:
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
			gd_table_extra = self.file_parser_gd_extra.run(gdpath)

		else:
			gd_df = None
			gd_table_extra = None

		if vcfpath:
			vcf_df = parse_vcf_file(vcfpath, set_index = self._set_table_index)
		else:
			vcf_df = None
		# Make sure the sequence id was identified properly.
		if index_df.index[0][0] == 'chrom1':
			# The sequence Id was not in the index file, but may be available in the gd or vcf files
			if gd_df is not None:
				available_table = gd_df
			elif vcf_df is not None:
				available_table = vcf_df
			else:
				available_table = None
			if available_table is not None:
				index_df = self.fix_index(index_df, available_table)
		# Merge the tables together.
		variant_df = self.merge_tables(index_df, gd_df, vcf_df)
		# Add the `sampleId` and `sampleName` columns
		variant_df[IsolateTableColumns.sample_id] = sample_id
		variant_df[IsolateTableColumns.sample_name] = sample_name
		return variant_df, coverage_df, junction_df, gd_table_extra
	def fix_index(self, source:pandas.DataFrame, key:pandas.DataFrame)->pandas.DataFrame:
		""" Fixes the index of `source` using the index from `key`"""

		if len(source) == len(key):
			# Can probably just use the indecies directly.
			# OK, so the index was renamed without throwing an exception.
			# Don't make a multiindex object. Reset the index, overwrite the original column values, then reindex.
			t = source.reset_index()
			key = key.reset_index()
			# Use tolist() to avoid errors if pandas decides to complain about the index.
			t[IsolateTableColumns.sequence_id] = key[IsolateTableColumns.sequence_id].tolist()
			t[IsolateTableColumns.position] = key[IsolateTableColumns.position].tolist()
			t = t.set_index([IsolateTableColumns.sequence_id, IsolateTableColumns.position])
		else:
			# No point in bothering.
			t = source
		return t


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
		# TODO: Make sure there are no duplicate column labels. For example, both the vcf and index tables have an 'alt' column
		# May need these if the gd or vcf files are missing.
		mutation_column = index[IsolateTableColumns.mutation].values
		annotation_column = index[IsolateTableColumns.annotation].values
		# For some reason the gd file appends 'chrom' to the sequence index.

		if gd is not None:
			# We don't care about most of the columns
			#logger.warning(list(gd.columns))
			reduced_gd = gd[self.gd_columns]
			variant_df: pandas.DataFrame = index.merge(reduced_gd, how = 'left', left_index = True, right_index = True)
		else:
			variant_df = index

			categories = [catagorize_mutation(a, m) for a, m in zip(mutation_column, annotation_column)]
			variant_df[IsolateTableColumns.mutation_category] = categories
			# Not generating the codon or amino acid columns since they aren't really used.
			variant_df['aminoAlt'] = variant_df['aminoRef'] = variant_df['codonAlt'] = variant_df['codonRef'] = 'unknown'

		if vcf is not None:
			# Remove the duplicate columns
			duplicate_columns = set(variant_df.columns) & set(vcf.columns)
			for extra_column in duplicate_columns:
				vcf.pop(extra_column)
			variant_df = variant_df.merge(vcf, how = 'left', left_index = True, right_index = True)
		else:
			if IsolateTableColumns.ref not in variant_df.columns:
				refs = [get_reference_from_mutation(m, a) for m, a in zip(mutation_column, annotation_column)]
				variant_df[IsolateTableColumns.ref] = refs
			if IsolateTableColumns.alt not in variant_df.columns:
				variant_df[IsolateTableColumns.alt] = variant_df[IsolateTableColumns.mutation].apply(get_alternate_from_mutation)
		return variant_df

	@staticmethod
	def get_summary(filename: Path, sample_id: str, sample_name: Optional[str] = None) -> Dict[str, Any]:
		return parse_summary_file(filename, sample_id, sample_name)

