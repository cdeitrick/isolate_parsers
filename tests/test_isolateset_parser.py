from pathlib import Path

from isolateset_parser import IsolateSetWorkflow
import pytest
data_folder = Path(__file__).with_name('data')


def test_parse_sample_map():
	sample_map_file = data_folder / "test_sample_map.txt"

	truth_data = {
		'sample1': 'sampleA',
		'sample2': 'sampleB',
		'sample3': 'sampleC',
		'sample4': '123'
	}
	test_data = IsolateSetWorkflow._parse_sample_map(sample_map_file)
	assert test_data == truth_data

@pytest.mark.parametrize(
	"value, expected",
	[
		("a,b,c", ["a", "b", "c"]),
		(None, []),
		("", []),
		("Tsv|456", ["Tsv|456"]),
		(["b", "d", "f"], ["b", "d", "f"])
	]
)
def test_parse_commandline_list(value, expected):
	result = IsolateSetWorkflow._parse_commandline_list(value)
	assert result == expected
