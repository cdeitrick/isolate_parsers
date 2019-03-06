import pytest
import pandas
import dataio
import datatools
@pytest.fixture
def variant_table()->pandas.DataFrame:
	string = """
		seq id	position	reference	E-01	E-02	E-04	E-05	E-06	E-07	E-08
		NODE_1	6397	G	A	A	G	A	A	A	A
		NODE_1	8463	G	G	G	G	G	G	G	G
		NODE_1	28954	A	G	G	A	G	G	G	G
		NODE_1	30052	A	G	G	A	G	G	G	G
		NODE_1	30921	C	C	C	C	C	C	C	A
		NODE_1	42817	G	G	G	G	G	G	G	G
		NODE_1	45057	C	T	T	T	T	T	T	T
		NODE_1	54213	C	T	T	T	T	T	T	T
		NODE_1	71269	G	A	A	G	A	A	A	G
		NODE_1	72044	C	C	C	C	C	C	C	C
		NODE_1	85276	G	G	G	G	G	G	G	G
		NODE_1	87791	C	T	T	C	T	T	T	T
		NODE_1	90358	C	T	T	C	T	T	T	T
		NODE_1	95837	C	T	T	C	T	T	T	T
		NODE_1	95879	G	G	G	G	G	G	G	G
		NODE_1	96380	G	G	G	G	G	G	G	G
		NODE_1	104532	C	T	T	C	T	T	T	T
		NODE_1	105935	C	C	C	C	C	C	T	C
		NODE_1	107432	T	T	T	T	T	T	T	C
		NODE_1	110679	C	T	T	T	T	T	T	T
	"""
	table = dataio.import_table(string)
	return table.set_index(['seq id', 'position'])

def test_filter_variants_in_all_samples(variant_table):
	expected_index = [
			("NODE_1",	6397),
			("NODE_1",	8463),
			("NODE_1",	28954),
			("NODE_1",	30052),
			("NODE_1",	30921),
			("NODE_1",	42817),
			("NODE_1",	71269),
			("NODE_1",	72044),
			("NODE_1",	85276),
			("NODE_1",	87791),
			("NODE_1",	90358),
			("NODE_1",	95837),
			("NODE_1",	95879),
			("NODE_1",	96380),
			("NODE_1",	104532),
			("NODE_1",	105935),
			("NODE_1",	107432)
	]
	result = datatools.filter_variants_in_all_samples(variant_table, 'reference')

	assert expected_index == list(result.index)

def test_variant_count(variant_table):
	expected = """
			seq id	position	count
		NODE_1	6397	6
		NODE_1	8463	0
		NODE_1	28954	6
		NODE_1	30052	6
		NODE_1	30921	1
		NODE_1	42817	0
		NODE_1	45057	7
		NODE_1	54213	7
		NODE_1	71269	5
		NODE_1	72044	0
		NODE_1	85276	0
		NODE_1	87791	6
		NODE_1	90358	6
		NODE_1	95837	6
		NODE_1	95879	0
		NODE_1	96380	0
		NODE_1	104532	6
		NODE_1	105935	1
		NODE_1	107432	1
		NODE_1	110679	7
	"""
	expected = dataio.import_table(expected)
	expected = expected.set_index(['seq id', 'position'])['count']
	reference = variant_table.pop('reference')
	result = datatools.variant_count(variant_table, reference)

	assert expected.to_dict() == result.to_dict()