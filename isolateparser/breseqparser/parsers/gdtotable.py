from pathlib import Path
from typing import *
import math
import pandas
from loguru import logger
from pprint import pprint


class GDToTable:
	def __init__(self):
		self.mutation_types = ["SNP", "SUB", "DEL", "INS", "MOB", "AMP", "CON", "INV", "AMP"]
		self.evidence_types = ["RA", "MC", "JC"]

		# There are a lot of columns to sort through.
		# Group columns together to make the table a little easier to read.
		self.field_groups = {
			'primary':   ['category_id', 'evidence_id', 'parent_ids'],
			'coverage':  [],
			'omit':      ['key', 'html_gene_name', 'locus_tags_inactivated'],
			'codon':     ['codon_new_seq', 'codon_number', 'codon_position', 'codon_ref_seq'],
			'amino':     ['aa_new_seq', 'aa_position', 'aa_ref_seq'],
			'gene':      [
				'locus_tag', 'locus_tags_overlapping', 'gene_name', 'gene_position', 'gene_product', 'gene_strand',
				'genes_overlapping',
				'genes_inactivated',
			],
			'mutation':  ['mutation_category', 'new_seq', 'seq_id', 'position', 'snp_type'],
			'insertion': ['repeat_length', 'repeat_new_copies', 'repeat_ref_copies', 'repeat_seq'],
			'deletion':  ['size']
		}

	@staticmethod
	def read_gd_file(filename: Path) -> List[str]:
		"""
		Reads in a gd file and returns non-commented lines.
		"""
		lines = filename.read_text().split('\n')
		lines = [i for i in lines if not i.startswith('#')]
		# remove empty lines
		lines = [i for i in lines if i.strip()]
		return lines

	def is_mutation(self, text: str) -> bool:
		""" Determines whether the text refers to a mutation or not."""
		key = text[:3]
		return key in self.mutation_types

	def is_evidence(self, text: str) -> bool:
		key = text[:2]
		return key in self.evidence_types

	@staticmethod
	def parse_keyword_fields(line: Union[str, List[str]]) -> Dict[str, str]:
		""" Splits the xpecified gd lines into a dictionary of key-value pairs. """
		if isinstance(line, str):
			line = line.split("\t")
		fields = [i.strip() for i in line]
		fields = [i.partition('=') for i in fields]
		# Remove the delimiter used to partition the field
		fields = [(i[0], i[2]) for i in fields]

		# The coverage values are presented as nn/nn, which excel interprets as a date.
		# Need to add quotes around those values so excel doesn't automatically convert them.

		# fields = [(i if '/' not in i[1] else (i[0], f"'{i[1]}'")) for i in fields]
		coverage_fields = [i for i in fields if ('cov' in i[0] and '/' in i[1])]
		other_fields = [i for i in fields if i not in coverage_fields]
		modified_coverage_fields = list()
		for key, value in coverage_fields:
			forward, reverse = value.split('/')
			numerator_key = f'{key}_forward'
			denominator_key = f"{key}_reverse"

			modified_coverage_fields.append((numerator_key, forward))
			modified_coverage_fields.append((denominator_key, reverse))

		fields = sorted(other_fields + modified_coverage_fields)
		return dict(fields)

	def parse_positional_fields(self, line: List[str]) -> Dict[str, str]:
		""" Parses the positional values of a line in the gd file."""
		line_type, evidence_id, parent_ids, seq_id, position = line[:5]
		data = {
			'category_id': line[0],
			'evidence_id': line[1],
			'parent_ids':  line[2],
		}
		# Mutations have more positional arguments than evidence fields, so `line` will have to be
		# modified accordingly.
		is_mutation = self.is_mutation(line_type)
		if is_mutation:
			data['seq_id'] = line[3]
			data['position'] = line[4]

		# Now update `data` the category-specific keyword arguments
		if line_type in {'SNP', 'INS'}:
			update = _apply_keywords_to_position(['new_seq'], [line[5]])
		elif line_type == 'DEL':
			update = _apply_keywords_to_position(['size'], [line[5]])
		elif line_type == 'AMP':
			update = _apply_keywords_to_position(['new_copy_number'], [line[5]])
		elif line_type == 'SUB':
			update = _apply_keywords_to_position(['size', 'new_seq'], [line[5], line[6]])
		elif line_type == 'MOB':
			update = _apply_keywords_to_position(['repeat_name', 'strand', 'duplication_size'], line[5:8])
		elif line_type == 'RA':
			keywords = ['seq_id', 'position', 'insert_position', 'ref_base', 'new_base']
			update = _apply_keywords_to_position(keywords, line[3:8])
		elif line_type == 'MC':
			update = _apply_keywords_to_position(['seq_id', 'start', 'end', 'start_range', 'end_range'], line[3:8])
		elif line_type == 'JC':
			keywords = [
				'side_1_seq_id', 'side_1_position', 'side_1_strand', 'side_2_seq_id',
				'side_2_position', 'side_2_strand', 'overlap'
			]
			update = _apply_keywords_to_position(keywords, line[3:10])
		elif line_type == 'UN':
			update = _apply_keywords_to_position(['seq_id', 'start', 'end'], line[3:6])
		else:
			message = f"Not a valid mutation or evidence type: '{line_type}'"
			raise ValueError(message)
		data.update(update)
		return data

	def parse_positional_fields2(self, line: List[str]) -> Dict[str, str]:
		""" Parses the positional values of a line in the gd file."""
		line_type, evidence_id, parent_ids, seq_id, position = line[:5]
		data = {
			'category_id': line[0],
			'evidence_id': line[1],
			'parent_ids':  line[2],
		}
		if self.is_mutation(line_type):
			data['seq_id'] = line[3]
			data['position'] = line[4]

		if line_type == 'SNP':
			data['new_seq'] = line[5]
		elif line_type == 'INS':
			data['new_seq'] = line[5]
		elif line_type == 'DEL':
			data['size'] = line[5]
		elif line_type == 'SUB':
			data['size'] = line[5]
			data['new_seq'] = line[6]
		elif line_type == 'MOB':
			data['repeat_name'] = line[5]
			data['strand'] = line[6]
			data['duplication_size'] = line[7]
		elif line_type == 'RA':
			data['seq_id'] = line[3]
			data['position'] = line[4]
			data['insert_position'] = line[5]
			data['ref_base'] = line[6]
			data['new_base'] = line[7]
		elif line_type == 'MC':
			data['seq_id'] = line[3]
			data['start'] = line[4]
			data['end'] = line[5]
			data['start_range'] = line[6]
			data['end_range'] = line[7]
		elif line_type == 'JC':
			data['side_1_seq_id'] = line[3]
			data['side_1_position'] = line[4]
			data['side_1_strand'] = line[5]
			data['side_2_seq_id'] = line[6]
			data['side_2_position'] = line[7]
			data['side_2_strand'] = line[8]
			data['overlap'] = line[9]
		elif line_type == 'AMP':
			data['new_copy_number'] = line[3]
		elif line_type == 'UN':
			data['seq_id'] = line[3]
			data['start'] = line[4]
			data['end'] = line[5]
		else:
			message = f"Not a valid mutation or evidence type: '{line_type}'"
			raise ValueError(message)

		return data

	def parse_line(self, line: str) -> Dict[str, str]:
		line = line.split("\t")

		positional_fields = self.parse_positional_fields(line)

		keyword_fields = line[len(positional_fields):]

		keyword_fields = self.parse_keyword_fields(keyword_fields)

		final_data = {**positional_fields, **keyword_fields}
		_line_type = line[0]
		final_data['category_id'] = _line_type

		return final_data

	def sort_columns(self, table: pandas.DataFrame) -> pandas.DataFrame:
		""" Groups columns together to make the table a little more legible."""
		groups = list()
		group_keys = ["primary", "mutation", 'deletion', "gene", "coverage", 'codon', 'amino', 'insertion']
		for key in group_keys:
			column_group = self.field_groups[key]
			group = [i for i in column_group if i in table.columns]
			groups += sorted(group)
		# Filter out useless columns.
		groups = [i for i in groups if i not in self.field_groups['omit']]
		other_columns = [i for i in table.columns if i not in groups]
		table = table[groups + other_columns]
		return table

	def run(self, filename: Path):
		lines = self.read_gd_file(filename)

		mutations = [i for i in lines if self.is_mutation(i)]
		evidence = [i for i in lines if self.is_evidence(i)]

		mutations = [self.parse_line(i) for i in mutations]
		evidence = [self.parse_line(i) for i in evidence]

		# Link evidence to mutation

		table_mutations = pandas.DataFrame(mutations)
		table_mutations.pop('evidence_id')
		table_evidence = pandas.DataFrame(evidence)
		table_evidence = table_evidence[[i for i in table_evidence.columns if i not in table_mutations.columns]]
		# Some mutations have multiple evidences. Split these off for now.
		table_mutations_single = table_mutations[~table_mutations['parent_ids'].str.contains(',')]

		merged_table = table_mutations_single.merge(
			table_evidence, left_on = 'parent_ids', right_on = 'evidence_id', how = 'left'
		)
		merged_table = merged_table.sort_values(by = ["mutation_category", 'seq_id', 'position'])
		merged_table = self.sort_columns(merged_table)
		return merged_table


def _apply_keywords_to_position(fields: List[str], values: List[Any]) -> Dict[str, Any]:
	row = {k: v for k, v in zip(fields, values)}
	return row


def main():
	filename = Path(
		"/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/lipuma/pipelines/SC1360/AU1064/breseq/output/evidence/annotated.gd")

	parser = GDToTable()

	result = parser.run(filename)

	result.info()
	result.to_csv("/home/cld100/Documents/sandbox/testmutationtable.tsv", sep = "\t", index = False)


if __name__ == "__main__":
	main()
