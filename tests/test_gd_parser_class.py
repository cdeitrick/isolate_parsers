from pathlib import Path

import pytest

from isolateparser.resultparser.parsers import parse_gd


@pytest.fixture
def gdparser() -> parse_gd.GDParser:
	return parse_gd.GDParser()


@pytest.fixture
def annotated_filename() -> Path:
	source = Path(__file__).parent / "data" / "gdfiles" / "annotated.gd"
	return source


@pytest.fixture
def annotated_gd_row() -> str:
	string = "SNP	3	44	2	1595623	A	aa_new_seq=I	aa_position=112	aa_ref_seq=F	" \
			 "codon_new_seq=ATC	codon_number=112	codon_position=1	codon_ref_seq=TTC	" \
			 "[db_xref=InterPro:IPR003754,InterPro:IPR007470]	" \
			 "[protein=protein of unknown function DUF513, hemX]	" \
			 "[protein_id=ABK09174.1][location=complement(2688435..26904 05)]	[gbkey=CDS]	" \
			 "snp_type=nonsynonymous"
	return string


def test_parse_row(gdparser, annotated_gd_row):
	expected = {
		"rowType":        "SNP",
		"rowId":          "3",
		"parentIds":      "44",
		"seq_id":          "2",
		"position":       "1595623",
		"new_seq":        "A",
		"aa_new_seq":     "I",
		"aa_position":    "112",
		"aa_ref_seq":     "F",
		"db_xref":        "InterPro:IPR003754,InterPro:IPR007470",
		'codon_new_seq':  'ATC',
		'codon_number':   '112',
		'codon_position': '1',
		'codon_ref_seq':  'TTC',
		"protein":        "protein of unknown function DUF513, hemX",
		"protein_id":     "ABK09174.1",
		"gbkey":          "CDS",
		"snp_type":       "nonsynonymous"
	}

	result = gdparser._parse_row(annotated_gd_row.split('\t'))

	assert result == expected


@pytest.mark.parametrize(
	"value,expected",
	[
		("codon_new_seq=ATC", "codon_new_seq=ATC"),
		("[protein_id=ABK09174.1][location=complement(2688435..26904 05)]", "protein_id=ABK09174.1"),
		("[gbkey=CDS]", "gbkey=CDS")
	]
)
def test_parse_keyword_argument(gdparser, value, expected):
	result = gdparser._parse_keyword_argument(value)

	assert result == expected


def test_convert_to_lines(gdparser):
	string = "abc\tdef\n123\nddd\tsss"
	expected = [["abc", "def"], ["123"], ["ddd", "sss"]]

	assert gdparser._convert_to_lines(string) == expected
