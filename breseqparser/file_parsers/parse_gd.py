from pathlib import Path
from typing import Dict, List, NamedTuple, Tuple, Union

import pandas
from dataclasses import dataclass, field


class _GDColumns(NamedTuple):
	description: str = 'description'
	gene: str = 'gene'
	mutation: str = 'mutation'
	position: str = 'position'
	sequence_id: str = 'seq id'
	alternate_base: str = 'baseAlt'
	#reference_base: str = 'baseRef'
	alternate_amino: str = 'aminoAlt'
	reference_amino: str = 'aminoRef'
	alternate_codon: str = 'codonAlt'
	reference_codon: str = 'codonRef'
	locus_tag: str = 'locusTag'
	mutation_category: str = 'mutationCategory'


GDColumns = _GDColumns()

MUTATION_KEYS: Dict[str, Dict[str, List[str]]] = {
	'snp': {
		'position': ['seqId', 'position', 'new_seq'],
		'keyword':  []
	},
	'sub': {
		'position': ['seqId', 'position', 'size', 'new_seq'],
		'keyword':  []
	},
	'del': {
		'position': ['seqId', 'position', 'size'],
		'keyword':  ['mediated', 'between', 'repeat_seq', 'repeat_length', 'repeat_ref_num',
			'repeat_new_copies']
	},
	'ins': {
		'position': ['seqId', 'position', 'new_seq'],
		'keyword':  ['repeat_seq', 'repeat_length', 'repeat_ref_num', 'repeat_new_copies', 'insert_position']
	},
	'mob': {
		'position': ['seqId', 'position', 'repeat_name', 'strand', 'duplication_size'],
		'keyword':  ['del_start', 'del_end', 'ins_start', 'ins_end', 'mob_region']
	},
	'amp': {
		'position': ['seqId', 'position', 'size', 'new_copy_number'],
		'keyword':  ['between', 'mediated', 'mob_region']
	},
	'con': {
		'position': ['seqId', 'position', 'size', 'region'],
		'keyword':  []
	},
	'inv': {
		'position': ['seqId', 'position', 'size'],
		'keyword':  []
	}
}
EVIDENCE_KEYS: Dict[str, Dict[str, List[str]]] = {
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


@dataclass
class Mutation:
	type: str
	id: str
	parentId: str
	seqId: str
	position: int

	details: Dict[str, Union[str, int]] = field(compare = False)

	def get(self, item: str, default = None) -> Union[str, int]:
		return self.details.get(item, default)

	@property
	def key(self) -> Tuple[str, int]:
		return self.seqId, int(self.position)


@dataclass
class Evidence:
	type: str
	id: str
	parentId: str
	details: Dict[str, Union[str, int]] = field(repr = False, compare = False)

	def get(self, item: str, default = None) -> Union[str, int]:
		return self.details.get(item, default)

	@property
	def seqId(self) -> str:
		result = self.details.get('seq_id', self.details.get('side_1'))
		return result

	@property
	def position(self) -> int:
		result = self.details.get('position')
		if result:
			result = int(result)
		return result

	@property
	def key(self) -> Tuple[str, int]:
		return self.seqId, self.position


def _sort_gd_file_rows(io: Union[str, Path]) -> Tuple[List[List[str]], List[List[str]]]:
	""" Reads a gd file and separates it into mutation and evidence tables.
		Parameters
		----------
		io: Union[str, Path]
			Either the text of a gd file or a path to a gd file.

		Returns
		-------
		mutations, evidence
	"""
	if isinstance(io, str):
		reader = io.split('\n')
	else:
		reader = io.read_text().split("\n")
	lines = [line.split('\t') for line in reader]
	mutation = filter(lambda s: s[0].lower() in MUTATION_KEYS, lines)
	evidence = filter(lambda s: s[0].lower() in EVIDENCE_KEYS, lines)
	return list(mutation), list(evidence)


def _get_row_position_and_sequence(row_type: str, other: List[str]) -> Dict[str, str]:
	"""
		Retrieves the `seqId` and `position` values from a gd file row, as well as basic information
		about the mutational change. For snps, this will just be the new base character, for
		a substitution, this will be the length and replacement sequence, etc.
	Parameters
	----------
	row_type: str
	other: List[str]
		The remainder of the gd file row after dropping the mutation type, id, and parent id.

	Returns
	-------
	Dict[str,str]
		A dictionary mapping the position and new sequence information to keys.
	"""
	if row_type in MUTATION_KEYS:
		mutation_config = MUTATION_KEYS[row_type]
	else:
		mutation_config = EVIDENCE_KEYS[row_type]

	position_keys = mutation_config['position']
	position_values = other[:len(position_keys)]
	position_map = dict(zip(position_keys, position_values))
	return position_map


def parse_annotated_gd_file_row(row: List[str]) -> Union[Mutation, Evidence]:
	""" Converts arow from a gd file into a `Mutation` or `Evidence` object."""
	row_type, row_id, parent_ids, *other = row
	row_type = row_type.lower()

	position_map = _get_row_position_and_sequence(row_type, other)
	keyword_values = dict([i.split('=') for i in other if '=' in i])

	if row_type in MUTATION_KEYS:
		r = Mutation(
			row_type, row_id, parent_ids,
			position_map['seqId'], int(position_map['position']),  # seq_id and position
			{**position_map, **keyword_values}
		)
	else:
		r = Evidence(
			type = row_type,
			id = row_id,
			parentId = parent_ids,
			details = {**position_map, **keyword_values}
		)
	return r


def parse_annotated_gd_file(path: Path) -> List[Union[Mutation, Evidence]]:
	""" Extracts information from the annotated gd file."""
	gd_table_mutations, gd_table_evidence = _sort_gd_file_rows(path)
	gd_table = gd_table_mutations + gd_table_evidence
	gd_data = map(parse_annotated_gd_file_row, gd_table)
	return list(gd_data)


def _extract_reference_base_from_codon(codon: str, position: str) -> str:
	try:
		reference_base = codon[int(position) - 1]
	except (ValueError, TypeError):
		reference_base = ""
	return reference_base


def generate_mutation_table(mutations: List[Mutation]) -> pandas.DataFrame:
	table = list()
	for mutation in mutations:
		position = int(mutation.get("position"))

		codon_ref = mutation.get('codon_ref_seq')
		codon_position = mutation.get('codon_position')
		# The reference base should be located from the vcf file since the gd file only has the reference codon.
		#reference_base = _extract_reference_base_from_codon(codon_ref, codon_position)

		row = {
			GDColumns.description:       mutation.get('gene_product'),
			GDColumns.gene:              mutation.get("gene_name"),
			GDColumns.mutation:          '',
			GDColumns.position:          position,
			GDColumns.sequence_id:       mutation.seqId,
			GDColumns.alternate_base:    mutation.get('new_seq'),
			#GDColumns.reference_base:    reference_base,
			GDColumns.alternate_amino:   mutation.get('aa_new_seq'),
			GDColumns.reference_amino:   mutation.get('aa_ref_seq'),
			GDColumns.locus_tag:         mutation.get('locus_tag'),
			GDColumns.mutation_category: mutation.get('mutation_category'),
			GDColumns.alternate_codon:   mutation.get('codon_new_seq'),
			GDColumns.reference_codon:   codon_ref
		}
		table.append(row)

	df = pandas.DataFrame(table)
	df = df[list(GDColumns)]
	# df.columns = GDColumns
	return df


def _extract_mutations(data: List[Union[Evidence, Mutation]]) -> List[Mutation]:
	return [i for i in data if isinstance(i, Mutation)]


def get_gd_filename(path: Path) -> Path:
	if path.is_dir():
		result = path / "output" / "evidence" / "annotated.vcf"
		if not result.exists():
			candidates = list(path.glob("**/annotated.gd"))
			if len(candidates) != 1:
				message = f"Invalid vcf file path: {path}"
				raise FileNotFoundError(message)
			result = candidates[0]
	elif path.suffix == '.gd':
		result = path
	else:
		message = f"Invalid GD path: {path}"
		raise ValueError(message)
	return result


def parse_gd_file(path: Path, set_index: bool = True) -> pandas.DataFrame:
	filename = get_gd_filename(path)
	gd_data = parse_annotated_gd_file(filename)
	mutations = _extract_mutations(gd_data)
	gd_df = generate_mutation_table(mutations)
	if set_index:
		gd_df.set_index(keys = [GDColumns.sequence_id, GDColumns.position], inplace = True)
	return gd_df


if __name__ == "__main__":
	_gd_file = Path(__file__).parent.parent.parent / "data" / "breseq_run" / "AU0074" / "breseq_output" / "output" / "evidence" / "annotated.gd"

	_gd_table = parse_gd_file(_gd_file.absolute())

	print(_gd_table.to_string())
