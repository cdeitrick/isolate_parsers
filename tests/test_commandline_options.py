from pathlib import Path

import pandas
import pytest
from loguru import logger
from isolateset_parser import IsolateSetWorkflow

"""
	program_options = load_program_options()
	isolateset_workflow = IsolateSetWorkflow(
		whitelist = program_options.whitelist,
		blacklist = program_options.blacklist,
		sample_map = program_options.sample_map,
		sample_regex = program_options.regex,
		use_filter = program_options.use_filter,
		snp_categories = program_options.snp_categories,
		generate_fasta = program_options.generate_fasta
	)
	isolateset_workflow.run(program_options.folder, program_options.reference_label)
"""


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
			"large deletion",
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
		'locusTag':            ["[pQBR0001]–[pQBR0478]", "PFLU0045/PFLU0046", "PFLU0872/PFLU0873", "PFLU3154/PFLU3155", "PFLU3571", "PFLU3777",
			"PFLU4443"],
		'mutationCategory':    ["large_deletion", "small_indel", "small_indel", "small_indel", "snp_nonsynonymous", "snp_nonsynonymous",
			"small_indel"],
		'position':            [1, 45881, 985333, 3447986, 3959631, 4173231, 4908233],
		'presentIn':           [3, 3, 3, 3, 2, 1, 1],
		'presentInAllSamples': [True, True, True, True, False, False, False],
		'ref':                 ['N/A', "G", "C", "", "G", "A", "G"],
		'seq id':              ["NC_009444", "NC_012660", "NC_012660", "NC_012660", "NC_012660", "NC_012660", "NC_012660"]
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

	# Compare using pieces of the table first to

	pandas.testing.assert_frame_equal(result_df, table)


if __name__ == "__main__":
	pass
