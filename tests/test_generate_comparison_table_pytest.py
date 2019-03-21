import pytest
import pandas
from file_generators import generate_comparison_table
from dataio import import_table
@pytest.fixture
def mutation_group()->pandas.DataFrame:
	string = """
	seq id	position	Sample	annotation	description	evidence	gene	mutation
	REL606	3762741	5K.index.html	K662I (AAA-ATA)	bifunctional (p)ppGpp synthetase II/ guanosine-3',5'-bis pyrophosphate 3'-pyrophosphohydrolase	RA	spoT -	A-T
	REL606	3762741	10K.index.html	K662I (AAA-ATA)	bifunctional (p)ppGpp synthetase II/ guanosine-3',5'-bis pyrophosphate 3'-pyrophosphohydrolase	RA	spoT -	A-T
	REL606	3762741	20K.index.html	K662I (AAA-ATA)	bifunctional (p)ppGpp synthetase II/ guanosine-3',5'-bis pyrophosphate 3'-pyrophosphohydrolase	RA	spoT -	A-T
	REL606	3762741	2K.index.html	K662I (AAA-ATA)	bifunctional (p)ppGpp synthetase II/ guanosine-3',5'-bis pyrophosphate 3'-pyrophosphohydrolase	RA	spoT -	A-T
	REL606	3762741	15K.index.html	K662I (AAA-ATA)	bifunctional (p)ppGpp synthetase II/ guanosine-3',5'-bis pyrophosphate 3'-pyrophosphohydrolase	RA	spoT -	A-T
	"""
	t = import_table(string, index = ["seq id", "position"])
	return t

def test_extract_annotation_from_group(mutation_group):
	result = generate_comparison_table._extract_annotation_from_group(mutation_group)
	expected = "K662I (AAA-ATA)"
	assert result == expected