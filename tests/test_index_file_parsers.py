from pathlib import Path

import pandas
import pytest

from isolateparser.breseqoutputparser.parsers import parse_index

data_folder = Path(__file__).parent / "data"
index_folder = data_folder / "index_files"

@pytest.fixture
def variant_table_parser()->parse_index.VariantTableParser():
	return parse_index.VariantTableParser()

@pytest.fixture
def index_parser()->parse_index.IndexParser():
	return parse_index.IndexParser()

@pytest.mark.parametrize(
	"filename",
	list(data_folder.joinpath("index_files").iterdir())
)
def test_parser_finishes(index_parser, filename):
	snp_table, cov_table, jun_table = index_parser.run("", filename)
	assert isinstance(snp_table, pandas.DataFrame)
	assert isinstance(cov_table, pandas.DataFrame)
	assert isinstance(jun_table, pandas.DataFrame)



@pytest.mark.parametrize(
	"number,expected",
	[
		("0.234", 0),
		("1,234.04", 1234),
		("12.3", 12)
	]
)
def test_to_integer(number, expected):
	result = parse_index.to_integer(number)
	assert expected == result


def test_extract_variant_table_headers(variant_table_parser):
	expected_clone = ["evidence", "seq\xa0id", "position", "mutation", "annotation", "gene", "description"]

	expected_population = ["evidence", "seq\xa0id", "position", "mutation", "freq", "annotation", "gene", "description"]

	text = """
		<!--Mutation Predictions -->
		<p>
		<!--Output Html_Mutation_Table_String-->
		<table border="0" cellspacing="1" cellpadding="3">
		<tr><th colspan="7" align="left" class="mutation_header_row">Predicted mutations</th></tr><tr>
		<th>evidence</th>
		<th>seq&nbsp;id</th>
		<th>position</th>
		<th>mutation</th>
		<th>annotation</th>
		<th>gene</th>
		<th width="100%">description</th>
		</tr>
		
		<!-- Item Lines -->
	"""
	text = "\n".join(i.strip() for i in text.split('\n'))
	headers = variant_table_parser._extract_table_headers(text)
	assert headers == expected_clone

@pytest.mark.parametrize(
	"value, expected",
	[
		("16123", 16123),
		(3.14159, 3),
		("161,171", 161171),
		("978,345:1", 978345)
	]
)
def test_clean_value(variant_table_parser, value, expected):
	result = variant_table_parser._clean_position(value)
	assert result == expected