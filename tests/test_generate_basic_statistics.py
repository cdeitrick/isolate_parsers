import pandas
import pytest

import dataio
from isolateparser.file_generators import generate_basic_statistics


@pytest.fixture
def comparison_table() -> pandas.DataFrame:
	string = """
		E-21	E-22	E-23	E-24	mutationCategory	position	presentIn	presentInAllSamples	ref	seq id
		C	C	C	C	snp_nonsynonymous	9463	17	0	G	NODE_10
		G	G	G	G	snp_nonsynonymous	33196	8	0	G	NODE_10
		T	T	T	G	snp_intergenic	36801	1	0	T	NODE_10
		G	G	GG	G	small_indel	43689	1	0	G	NODE_10
		A	A	A	A	snp_intergenic	85184	21	1	G	NODE_10
		T	T	T	T	snp_intergenic	85235	3	0	T	NODE_10
	"""

	return dataio.import_table(string)


def test_calculate_basic_statistics(comparison_table):
	expected = """
		sampleName	small_indel	snp_intergenic	snp_nonsynonymous	dN/dS
		E-21	0	1	1	0
		E-22	0	1	1	0
		E-23	1	1	1	0
		E-24	0	2	1	0
	"""
	expected = dataio.import_table(expected)
	result = generate_basic_statistics.calculate_basic_statistics(comparison_table)
	# assert list(expected.columns) == list(result.columns)
	pandas.testing.assert_frame_equal(expected, result, check_dtype = False)
