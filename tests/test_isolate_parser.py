from pathlib import Path

import pytest
from loguru import logger

from isolateparser.breseqoutputparser import breseq_output_parser

data_folder = Path(__file__).parent / 'data'


@pytest.fixture
def breseq_clone() -> Path:
	# In case I need to customize the Path object
	clone_folder = data_folder / "Clonal_Output" / "breseq_output"
	return clone_folder


@pytest.fixture
def breseq_isolate_parser() -> breseq_output_parser.BreseqOutputParser:
	return breseq_output_parser.BreseqOutputParser(use_filter = False)


def test_get_sample_name():
	result = breseq_output_parser.get_sample_name(Path("Clonal_Output"))
	assert result == 'Clonal_Output'

	result = breseq_output_parser.get_sample_name(Path("Clonal_Output") / "breseq output")
	assert result == 'Clonal_Output'


def test_parse_breseq_isolate_with_only_index_path(breseq_clone, breseq_isolate_parser):
	# TODO: refactor this so it doesn't read from IsolateTableColumns
	# TODO: test different combinations of index,vcf,gd paths.
	variant_table, coverage_table, junction_table = breseq_isolate_parser.run('testIsolate', breseq_clone / "output" / "index.html")
	logger.info(list(variant_table.columns))
	logger.info(list(breseq_output_parser.IsolateTableColumns))

	# Since the workflow was only provided the index file, it is missing the additional columns added from the vcf and gd files.

	# TODO: This test currently fails because the 'alt', 'ref', 'locustag', and 'mutationCategory' columns are missing due to the
	# missing vcf and gd files. Refactor the code so that these columns are parsed form the index file instead.
	assert list(variant_table.columns) == list(i for i in breseq_output_parser.IsolateTableColumns if i not in {'seq id', 'position'})
