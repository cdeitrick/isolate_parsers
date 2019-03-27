import pytest
from pathlib import Path

from breseqparser import isolate_parser

data_folder = Path(__file__).parent / 'data'

@pytest.fixture
def breseq_clone()->Path:
	# In case I need to customize the Path object
	clone_folder = data_folder / "Clonal_Output" / "breseq_output"
	return clone_folder


def test_get_sample_name():
	result = isolate_parser.get_sample_name(Path("Clonal_Output"))
	assert result == 'Clonal_Output'

	result = isolate_parser.get_sample_name(Path("Clonal_Output") / "breseq output")
	assert result == 'Clonal_Output'

def test_get_file_ocations(breseq_clone):
	expected_index = breseq_clone / "output" / "index.html"
	expected_vcf = breseq_clone / "data" / "output.vcf"
	expected_gd = breseq_clone / "output" / "evidence" / "annotated.gd"

	index, vcf, gd = isolate_parser._get_file_locations(breseq_clone)

	assert index == expected_index
	assert vcf == expected_vcf
	assert gd == expected_gd

def test_parse_breseq_isolate(breseq_clone):
	variant_table, coverage_table, junction_table = isolate_parser.parse_breseq_isolate(breseq_clone, isolate_id = 'testIsolate')
	assert list(variant_table.columns) == list(isolate_parser.IsolateTableColumns)

