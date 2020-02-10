"""
	Implements the classes which indicate which columns are available from each table.
	These should be consolidated into a single file instead of beings spread theoughout the codebase.
"""
from typing import NamedTuple
class IsolateTableColumns(NamedTuple):
	# Defines the column names for the isolate table.
	# Mainly used as a reminder of what the final column labels should be.
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
IsolateTableColumns = IsolateTableColumns()

import enum

class Columns(enum.Enum):
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

if __name__ == "__main__":
	print(list(Columns))