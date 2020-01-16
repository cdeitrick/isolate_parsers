from pathlib import Path

import pytest
from loguru import logger

from isolateparser import breseqparser

data_folder = Path(__file__).parent / 'data'


@pytest.fixture
def breseq_clone() -> Path:
	# In case I need to customize the Path object
	clone_folder = data_folder / "Clonal_Output"
	return clone_folder


@pytest.fixture
def breseq_isolate_parser() -> breseqparser.BreseqFolderParser:
	return breseqparser.BreseqFolderParser(use_filter = False)


def test_get_sample_name():
	result = breseqparser.get_sample_name(Path("Clonal_Output"))
	assert result == 'Clonal_Output'

	result = breseqparser.get_sample_name(Path("Clonal_Output") / "breseq output")
	assert result == 'Clonal_Output'


def test_parse_breseq_isolate_with_only_index_path(breseq_clone, breseq_isolate_parser):
	# TODO: refactor this so it doesn't read from IsolateTableColumns
	# TODO: test different combinations of index,vcf,gd paths.
	variant_table, coverage_table, junction_table = breseq_isolate_parser.run('testIsolate', breseq_clone / "output" / "index.html")
	logger.info(list(variant_table.columns))
	logger.info(list(breseqparser.IsolateTableColumns))

	# Since the workflow was only provided the index file, it is missing the additional columns added from the vcf and gd files.

	# TODO: This test currently fails because the 'alt', 'ref', 'locustag', and 'mutationCategory' columns are missing due to the
	# missing vcf and gd files. Refactor the code so that these columns are parsed form the index file instead.
	expected = ['sampleId', 'sampleName', 'annotation', 'description', 'frequency', 'gene', 'mutationCategory', 'alt', 'ref', 'mutation']
	extra = ['aminoRef', 'aminoAlt', 'codonRef', 'codonAlt']
	assert sorted(variant_table.columns) == sorted(expected + extra)

@pytest.mark.parametrize(
	"value, annotation, expected",
	[
		("D425,094 bp", "", 				"."),
		("+G", 			"intergenic (+61/+24)",		"."),
		("+C", 			"intergenic (+16/‑137)", 	"."),
		("C-T",			"intergenic (+112/+695)", 	"C"),
		("G-C", 		"intergenic (+552/+255)", 	"G"),
		("C-A", 		"S211A (TCG→GCG)", 			"C")
	]
)
def test_get_reference_from_mutation(value, annotation, expected):

	result = breseqparser.breseq_folder_parser.get_reference_from_mutation(value, annotation)

	assert result == expected

@pytest.mark.parametrize(
	"value, expected",
	[
		("D425,094 bp", "D425,094 bp"),
		("+G", "G"),
		("+C", "C"),
		("C-T", "T"),
		("G-C", "C"),
		("C-A", "A")
	]
)
def test_get_alternate_from_mutation(value, expected):

	result = breseqparser.breseq_folder_parser.get_alternate_from_mutation(value)

	assert result == expected