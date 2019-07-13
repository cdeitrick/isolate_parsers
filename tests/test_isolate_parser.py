from pathlib import Path

import pytest
from loguru import logger

from isolateparser.breseqparser import isolate_parser

data_folder = Path(__file__).parent / 'data'


@pytest.fixture
def breseq_clone() -> Path:
	# In case I need to customize the Path object
	clone_folder = data_folder / "Clonal_Output" / "breseq_output"
	return clone_folder


@pytest.fixture
def breseq_isolate_parser() -> isolate_parser.BreseqOutputParser:
	return isolate_parser.BreseqOutputParser(use_filter = False)


def test_get_sample_name():
	result = isolate_parser.get_sample_name(Path("Clonal_Output"))
	assert result == 'Clonal_Output'

	result = isolate_parser.get_sample_name(Path("Clonal_Output") / "breseq output")
	assert result == 'Clonal_Output'


def test_get_file_locations(breseq_clone):
	expected_index = breseq_clone / "output" / "index.html"
	expected_vcf = breseq_clone / "data" / "output.vcf"
	expected_gd = breseq_clone / "output" / "evidence" / "annotated.gd"

	index, vcf, gd = isolate_parser.BreseqOutputParser.get_file_locations(breseq_clone)

	assert index == expected_index
	assert vcf == expected_vcf
	assert gd == expected_gd


def test_parse_breseq_isolate(breseq_clone, breseq_isolate_parser):
	variant_table, coverage_table, junction_table = breseq_isolate_parser.run(breseq_clone, sample_id = 'testIsolate')
	logger.info(list(variant_table.columns))
	logger.info(list(isolate_parser.IsolateTableColumns))
	assert list(variant_table.columns) == list(i for i in isolate_parser.IsolateTableColumns if i not in {'seq id', 'position'})
