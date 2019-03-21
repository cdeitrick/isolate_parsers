import pytest
from pathlib import Path
from breseqparser.file_parsers import parse_index
import pandas

data_folder = Path(__file__).parent / "data"

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

def test_parse_snp_table():
	pass