from pathlib import Path

import pandas
import pytest

from isolateparser.breseqoutputparser.file_parsers import parse_index

data_folder = Path(__file__).parent / "data"
index_folder = data_folder / "index_files"



@pytest.mark.parametrize(
	"filename",
	list(data_folder.joinpath("index_files").iterdir())
)
def test_parser_finishes(filename):
	snp_table, cov_table, jun_table = parse_index.parse_index_file("", filename)
	assert isinstance(snp_table, pandas.DataFrame)
	assert isinstance(cov_table, pandas.DataFrame)
	assert isinstance(jun_table, pandas.DataFrame)


@pytest.mark.parametrize(
	"folder,expected",
	[
		(data_folder / "index_files" / "index.html", data_folder / "index_files" / "index.html"),
		(data_folder / "index_files" / "index (1).html", data_folder / "index_files" / "index (1).html")
	]
)
def test_get_index_filename(folder, expected):
	result = parse_index.get_index_filename(folder)
	assert result == expected


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


def test_extract_variant_table_headers():
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
	headers = parse_index._extract_variant_table_headers(text)
	assert headers == expected_clone
