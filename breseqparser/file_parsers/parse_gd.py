from pathlib import Path
from typing import List, Dict, Union, Tuple
import csv
from dataclasses import dataclass, field
import pandas
MUTATION_KEYS = {
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
EVIDENCE_KEYS = {
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


class GenomeDiffParser:
	"""
		Combines information from each data source describing the variants found in the isolate.
	Parameters
	----------
	path: Path
		path to the breseq output folder.

	"""

	def __init__(self, path: Union[str, Path], sample_id: str = None):
		if sample_id:
			self.sample_id = sample_id
		else:
			self.sample_id = [i for i in path.parts if i not in {'output', 'annotated', 'evidence'} and not i.endswith('.gd')][-1]
		path = Path(path)
		self.mutation_keys = MUTATION_KEYS
		self.evidence_keys = EVIDENCE_KEYS

		# The annotated gd_file contains essentially the same information as the index.html file.
		if path.is_dir():
			gd_file = self._search_for_annotated_gd_file(path)
		else:
			gd_file = path

		self.gd_data = self.parse_annotated_gd(gd_file)

	@staticmethod
	def _search_for_annotated_gd_file(path: Path) -> Path:
		basic = path / "output" / "evidence" / "annotated.gd"
		if not basic.exists():
			candidates = list(path.glob("**/annotated.gd"))
			if len(candidates) > 0:
				basic = candidates[0]
			else:
				message = "Cannot find the annotated gd file in the folder {}".format(path)
				raise FileNotFoundError(message)

		return basic

	@staticmethod
	def save(table: List[Dict], path: Path) -> None:
		import itertools
		fieldnames = itertools.chain.from_iterable([i.keys() for i in table])
		fieldnames = set(fieldnames)
		with path.open('w') as file1:
			writer = csv.DictWriter(file1, fieldnames = fieldnames, delimiter = "\t")
			writer.writeheader()
			for row in table:
				writer.writerow(row)

	def parse_annotated_gd(self, path: Path) -> List[Union[Mutation, Evidence]]:
		""" Extracts information from the annotated gd file."""

		gd_table_mutations, gd_table_evidence = self._separated_gd_rows(path)
		gd_table = gd_table_mutations + gd_table_evidence
		gd_data = list()
		for row in gd_table:
			row_type, row_id, parent_ids, *other = row
			row_type = row_type.lower()

			if row_type in self.mutation_keys:
				mutation_config = self.mutation_keys[row_type]
			else:
				mutation_config = self.evidence_keys[row_type]

			position_keys = mutation_config['position']

			position_values = other[:len(position_keys)]
			position_map = dict(zip(position_keys, position_values))
			keyword_values = dict([i.split('=') for i in other if '=' in i])

			if row_type in self.mutation_keys:
				r = Mutation(
					row_type, row_id, parent_ids,
					position_values[0], int(position_values[1]),  # seq_id and position
					{**position_map, **keyword_values}
				)
			else:
				r = Evidence(
					type = row_type,
					id = row_id,
					parentId = parent_ids,
					details = {**position_map, **keyword_values}
				)
			gd_data.append(r)
		return gd_data

	def _separated_gd_rows(self, filename: Path) -> Tuple[List[List[str]], List[List[str]]]:
		""" Separates a gd file into mutation and evidence tables."""
		mutations = list()
		evidence = list()
		with filename.open() as gd_file:
			reader = gd_file.read().split("\n")
			for line in reader:
				row = line.split("\t")
				first_element = row[0].lower()
				if first_element in self.mutation_keys.keys():
					mutations.append(row)
				elif first_element in self.evidence_keys.keys():
					evidence.append(row)
				elif first_element in {}:
					pass
				else:
					pass
		return mutations, evidence

	def mutations(self, kind = 'all') -> List[Mutation]:
		ls = list()
		for i in self.gd_data:
			if isinstance(i, Mutation) and (kind == 'all' or i.type == kind):
				ls.append(i)
		return ls

	def evidence(self, kind = 'all') -> List[Evidence]:
		ls = list()
		for i in self.gd_data:
			if isinstance(i, Evidence) and (kind == 'all' or i.type == kind):
				ls.append(i)
		return ls

	def generate_mutation_table(self)->pandas.DataFrame:
		table = list()
		for mutation in self.mutations():
			description = mutation.get('gene_product')
			gene_name = mutation.get("gene_name")
			position = int(mutation.get("position"))

			codon_ref = mutation.get('codon_ref_seq')
			codon_position = mutation.get('codon_position')
			try:
				reference_base = codon_ref[int(codon_position)-1]
			except (ValueError, TypeError):
				reference_base = ""

			row = {
				'sampleName': self.sample_id,
				'description': description,
				'gene': gene_name,
				'mutation': '',
				'position': position,
				'seq id': mutation.seqId,
				'baseAlt': mutation.get('new_seq'),
				'baseRef': reference_base,
				'aminoAlt': mutation.get('aa_new_seq'),
				'aminoRef': mutation.get('aa_ref_seq'),
				'locusTag': mutation.get('locus_tag'),
				'mutationCategory': mutation.get('mutation_category')
			}
			table.append(row)

		df = pandas.DataFrame(table)
		return df

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

def parse_gd(path:Path, sample_id:str)->pandas.DataFrame:
	filename = get_gd_filename(path)
	gd_data = GenomeDiffParser(filename, sample_id)
	return gd_data.generate_mutation_table()

if __name__ == "__main__":
	gd_file = Path(__file__).parent.parent.parent / "data" / "breseq_run" / "AU0074" / "breseq_output" / "output" / "evidence" / "annotated.gd"

	gd_table = parse_gd(gd_file.absolute(), 'AU0074')

	#print(gd_table.to_string())
