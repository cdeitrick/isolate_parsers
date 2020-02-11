from pathlib import Path
from typing import Dict, List, Tuple, Union

import pandas
from loguru import logger

MUTATION_FIELDS: Dict[str, Dict[str, List[str]]] = {
	'snp': {
		'position': ['seq_id', 'position', 'new_seq'],
		'keyword':  []
	},
	'sub': {
		'position': ['seq_id', 'position', 'size', 'new_seq'],
		'keyword':  []
	},
	'del': {
		'position': ['seq_id', 'position', 'size'],
		'keyword':  ['mediated', 'between', 'repeat_seq', 'repeat_length', 'repeat_ref_num',
			'repeat_new_copies']
	},
	'ins': {
		'position': ['seq_id', 'position', 'new_seq'],
		'keyword':  ['repeat_seq', 'repeat_length', 'repeat_ref_num', 'repeat_new_copies', 'insert_position']
	},
	'mob': {
		'position': ['seq_id', 'position', 'repeat_name', 'strand', 'duplication_size'],
		'keyword':  ['del_start', 'del_end', 'ins_start', 'ins_end', 'mob_region']
	},
	'amp': {
		'position': ['seq_id', 'position', 'size', 'new_copy_number'],
		'keyword':  ['between', 'mediated', 'mob_region']
	},
	'con': {
		'position': ['seq_id', 'position', 'size', 'region'],
		'keyword':  []
	},
	'inv': {
		'position': ['seq_id', 'position', 'size'],
		'keyword':  []
	}
}

EVIDENCE_FIELDS: Dict[str, Dict[str, List[str]]] = {
	'ra': {
		'position': ['seq_id', 'position', 'insert_position', 'ref_base', 'new_base'],
		'keyword':  []
	},
	'mc': {
		'position': ['seq_id,', 'start', 'end', 'start_range', 'end_range'],
		'keyword':  []
	},
	'jc': {
		'position': ['side_1'],
		'keyword':  []
	}
}
GD_FIELDS = {**MUTATION_FIELDS, **EVIDENCE_FIELDS}


class GDParser:
	#  ['aminoAlt', 'aminoRef', 'codonAlt', 'codonRef', 'locusTag', 'mutationCategory']
	column_map = {
		'seq_id': 'seq id',
		'aa_new_seq': 'aminoAlt',
		'aa_ref_seq': 'aminoRef',
		'codon_new_seq': 'codonAlt',
		'codon_ref_seq': 'codonRef',
		'mutation_category': 'mutationCategory'
	}
	def __init__(self):
		pass
	def run(self, io: Union[str, Path], set_index:bool = True) -> pandas.DataFrame:
		content = self._load_content(io)
		lines = self._convert_to_lines(content)

		# Convert each line to a dictionary.
		gd_data: List[Dict[str, str]] = [self._parse_row(l) for l in lines]

		df = pandas.DataFrame(gd_data)
		df = self._rename_columns(df)

		# We don't care about the evidence rows.
		df = df[~df['rowType'].isin(['RA', 'MC', 'JC'])]

		# Make sure the position column is actually an int rather than a str
		df['position'] = df['position'].astype(int)

		# The `seq id` field may be formatted a little differently than in the index file.
		#df['seq id'] = df['seq id'].apply(lambda s: 'chrom'+str(s))
		if set_index:
			df.set_index(keys = ['seq id', 'position'], inplace = True)
		return df

	@staticmethod
	def _convert_to_lines(content: str) -> List[List[str]]:
		""" Splits a file into lines while omitting unnecessary lines."""
		lines = content.split('\n')
		# Remove empty lines and the header
		lines = [i for i in lines if (bool(i) and not i.startswith('#'))]

		# Tokenize the line
		tokens = [i.split('\t') for i in lines]

		clean_tokens = [i for i in tokens if i[0].lower() != 'un']

		return clean_tokens

	@staticmethod
	def _categorize_lines(lines: List[List[str]]) -> Tuple[List[List[str]], List[List[str]]]:
		""" Groups lines according to whether they represent a `mutation` or `evidence`."""
		mutations = list()
		evidence = list()

		for line in lines:
			line_type = line[0]

			# Mutations are three characters, evidence types are two characters
			if len(line_type) == 2:
				evidence.append(line)
			elif len(line_type) == 3:
				mutations.append(line)
			else:
				message = f"Cannot identify whether this line is a `mutation` or `evidene`: {line}"
				logger.warning(message)
		return mutations, evidence

	@staticmethod
	def _load_content(io: Union[str, Path]) -> str:
		if isinstance(io, Path):
			content = io.read_text()
		else:
			try:
				content = Path(io).read_text()
			except OSError:
				# `io` is too long, because it is actually the contents of the file.
				content = io
		return content

	def _parse_row(self, line: List[str]) -> Dict[str, str]:
		""" Converts a line in the gd file into a dictionary mapping parameters to their corresponding values."""
		row_type = line[0].lower()
		type_definition = GD_FIELDS[row_type]
		position_fields = ['rowType', 'rowId', 'parentIds'] + type_definition['position']

		position_values = dict(zip(position_fields, line))  # zip() should stop when last item in `posiiton_files` is reached.
		keyword_values = dict()
		for keyword_argument in [i for i in line if '=' in i]:
			arg = self._parse_keyword_argument(keyword_argument)
			key, _, value = arg.partition('=')  # Use `partition` in case the value has a `=` character.
			keyword_values[key] = value

		return {**position_values, **keyword_values}

	@staticmethod
	def _parse_keyword_argument(argument: str) -> str:
		if argument.startswith('['): argument = argument[1:]
		if argument.endswith(']'): argument = argument[:-1]

		if '][' in argument:
			argument = argument.split('][')[0]  # Ignore 'location'
		return argument

	def _rename_columns(self, df:pandas.DataFrame)->pandas.DataFrame:
		""" Renames som eo fthe columns form the gd file."""

		for old_col, new_col in self.column_map.items():
			if old_col not in df.columns:
				logger.warning(f"Cannot find '{old_col}'")
				continue
			df[new_col] = df[old_col]
			del df[old_col]
		return df