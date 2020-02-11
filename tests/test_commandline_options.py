from pathlib import Path

import pandas
import pytest
from loguru import logger
from isolateset_parser import IsolateSetWorkflow
import math

@pytest.fixture
def table() -> pandas.DataFrame:
	data = {
		'S1_58BA':             ['G', 'GG', 'CC', '', 'A', 'A', 'G'],
		'S2_58BA':             ['G', 'GG', 'CC', '', 'A', 'A', 'G'],
		'SC_58SM':             ['G', 'GG', 'CC', '', 'G', 'C', 'GG'],
		'annotation':          [
			'D425,094 bp', 'intergenic (+65/+20)', 'intergenic (+17/-136)',
			'intergenic (+57/+21)', 'M350I (ATG-ATA)', 'T238P (ACC-CCC)',
			'coding (322/1476 nt)'
		],
		'description':         [
			"large_deletion",
			"putative lipoprotein/putative hydrolase",
			"microcin-processing peptidase 1. Unknown type peptidase. MEROPS family U62/hypothetical protein",
			"hypothetical protein/putative helicase",
			"putative GGDEF domain signaling protein",
			"hybrid sensory histidine kinase in two-component regulatory system with UvrY",
			"putative two-component system response regulator nitrogen regulation protein NR(I)"
		],
		'gene':                [
			"parA-pQBR0478", "PFLU0045 - / - PFLU0046",
			"PFLU0872 - / - PFLU0873",
			"PFLU3154 - / - PFLU3155",
			"PFLU3571 -", "PFLU3777 -",
			"PFLU4443 -"
		],
		#'locusTag':            ["[pQBR0001]–[pQBR0478]", "PFLU0045/PFLU0046", "PFLU0872/PFLU0873", "PFLU3154/PFLU3155", "PFLU3571", "PFLU3777",
	#	"PFLU4443"],
		'locusTag': [math.nan]*7,
		'mutationCategory':    ["large_deletion", "small_indel", "small_indel", "small_indel", "snp_nonsynonymous", "snp_nonsynonymous",
			"small_indel"],
		'position':            [1, 45881, 985333, 3447986, 3959631, 4173231, 4908233],
		'presentIn':           [3, 3, 3, 3, 2, 1, 1],
		'presentInAllSamples': [True, True, True, True, False, False, False],
		'ref':                 ['N/A', "G", "C", "", "G", "A", "G"],
		'seq id':              ["NC_009444", "NC_012660", "NC_012660", "NC_012660", "NC_012660", "NC_012660", "NC_012660"],
		"inReference":			[True,True,True,True, True,False,False]
	}

	df = pandas.DataFrame(data)

	return df


@pytest.fixture
def run_folder() -> Path:
	folder_data = Path(__file__).parent / "data"
	folder = folder_data / "set_output"
	return folder


def read_table(filename: Path) -> pandas.DataFrame:
	df = pandas.read_excel(filename, sheet_name = "variant comparison")
	return df


@pytest.fixture
def sample_group():
	pass


def test_default_commandline_arguments(run_folder, table):
	isolate_set_workflow = IsolateSetWorkflow()
	reference_label = "S1_58BA"
	result = isolate_set_workflow.run(run_folder, reference_label)

	result_df = read_table(result)

	# Compare using pieces of the table first since pandas.testing.assert_frame_equal doesn't report much.
	assert sorted(result_df.columns) == sorted(table.columns)
	assert sorted(result_df.index) == sorted(table.index)

	# Now test whether the column values are the same. Again, the pandas testing functions aren't helpful when doing this manually.
	for column in table.columns:
		if column in ["S1_58BA", "S2_58BA", "SC_58SM"]: continue
		logger.debug(column)
		assert result_df[column].tolist() == table[column].tolist()

	#pandas.testing.assert_frame_equal(result_df, table)


if __name__ == "__main__":
	pass
