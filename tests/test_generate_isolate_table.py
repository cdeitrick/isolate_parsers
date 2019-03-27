import pytest
import pandas
import dataio
from file_generators import generate_isolate_table
@pytest.fixture
def isolate_table()->pandas.DataFrame:
	string = """
	E-21	E-22	E-23	E-24	annotation	description	gene	locusTag	mutationCategory	position	presentIn	presentInAllSamples	ref	seq id
	C	C	C	C	A312P (GCG-CCG)	BPFJKKCJ_05542 -		snp_nonsynonymous	9463	17	0	G	NODE_10
	G	G	G	G	G65S (GGT-AGT)	BPFJKKCJ_05564 -		snp_nonsynonymous	33196	8	0	G	NODE_10
	T	T	T	G	intergenic (-131/-525)	BPFJKKCJ_05568 - / - BPFJKKCJ_05569	–/–	snp_intergenic	36801	1	0	T	NODE_10
	G	G	GG	G	intergenic (-32/+95)	BPFJKKCJ_05577 - / - BPFJKKCJ_05578	–/–	small_indel	43689	1	0	G	NODE_10
	A	A	A	A	intergenic (+188/-383)	BPFJKKCJ_05614 - / - BPFJKKCJ_05615	–/–	snp_intergenic	85184	21	1	G	NODE_10
	T	T	T	T	intergenic (+239/-332)	BPFJKKCJ_05614 - / - BPFJKKCJ_05615	–/–	snp_intergenic	85235	3	0	T	NODE_10
	"""
	return dataio.import_table(string)
def test_save_tables(isolate_table, tmp_path):
	excel_file = tmp_path / "test.xlsx"
	generate_isolate_table.save_isolate_table({'test': isolate_table}, excel_file)
	assert excel_file.exists()
