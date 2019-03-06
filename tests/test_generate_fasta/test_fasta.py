import pytest
from file_generators import generate_fasta
import pandas
@pytest.fixture
def combined_table()->pandas.DataFrame:
	string = """
		seq id	position	reference	E-01	E-02	E-05	E-06	E-07	E-09	E-10
		NODE_1	6397	G	A	A	A	A	A	A	A
		NODE_1	8463	G	G	G	G	G	G	G	G
		NODE_1	28954	A	G	G	G	G	G	G	G
		NODE_1	30052	A	G	G	G	G	G	G	G
		NODE_1	42817	G	G	G	G	G	G	G	G
		NODE_1	45057	C	T	T	T	T	T	T	T
		NODE_1	54213	C	T	T	T	T	T	T	T
		NODE_1	71269	G	A	A	A	A	A	A	A
		NODE_1	72044	C	C	C	C	C	C	C	C
		NODE_1	85276	G	G	G	G	G	G	G	G
		NODE_1	87791	C	T	T	T	T	T	T	T
		NODE_1	90358	C	T	T	T	T	T	T	T
		NODE_1	95837	C	T	T	T	T	T	T	T
		NODE_1	95879	G	G	G	G	G	G	G	G
		NODE_1	96380	G	G	G	G	G	G	G	G
		NODE_1	104532	C	T	T	T	T	T	T	T
		NODE_1	105935	C	C	C	C	C	T	C	C
		NODE_1	110679	C	T	T	T	T	T	T	T
		NODE_1	136698	G	G	A	A	A	A	A	A
		NODE_1	176069	G	A	A	A	A	A	A	A
		NODE_1	206410	A	G	G	G	G	G	G	G
		NODE_1	247696	G	G	G	G	G	G	G	G
		NODE_1	256818	C	T	T	T	T	T	T	T
		NODE_1	274515	C	T	T	T	T	T	T	T
		NODE_1	282503	C	T	T	T	T	T	T	T
		NODE_1	300752	G	G	G	G	G	G	G	G
		NODE_1	305443	A	G	G	G	G	G	G	G
		NODE_1	324904	A	G	G	G	G	G	G	G
		NODE_1	367587	A	G	G	G	G	G	G	G
	"""

def test_convert_combined_table_to_aligned_table():
	pass

def test_generate_reference_sequence():
	pass