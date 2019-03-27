
from pathlib import Path

from isolateset_parser import _parse_commandline_list, _parse_sample_map

data_folder = Path(__file__).with_name('data')


def test_parse_sample_map():
	sample_map_file = data_folder / "test_sample_map.txt"

	truth_data = {
		'sample1': 'sampleA',
		'sample2': 'sampleB',
		'sample3': 'sampleC',
		'sample4': '123'
	}
	test_data = _parse_sample_map(sample_map_file)
	assert test_data == truth_data


def test_parse_commandline_list():
	teststring = "a,b,c"
	assert teststring.split(',') == _parse_commandline_list(teststring)
