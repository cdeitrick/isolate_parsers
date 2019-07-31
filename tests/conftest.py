import pytest
from pathlib import Path
def checkdir(path)->Path:
	path = Path(path)
	if not path.exists():
		path.mkdir()
	return path
@pytest.fixture
def breseq_pipeline_output_empty(tmp_path)->Path:
	""" Sets up multiple folders each with its own breseq output.

	"""

	parent_folder = checkdir(tmp_path / "parent_folder")

	sample_1_folder = checkdir(parent_folder / "sample1")
	sample_1_folder_output = checkdir(sample_1_folder / "output")
	sample_1_folder_data = checkdir(sample_1_folder / "data")

	sample_2_folder = checkdir(parent_folder / "sample2")
	sample_2_folder_breseq = checkdir(sample_2_folder / "breseq")
	sample_2_folder_output = checkdir(sample_2_folder_breseq / "output")
	sample_2_folder_data = checkdir(sample_2_folder_breseq / "data")

	sample_3_folder = checkdir(parent_folder / "AU1234_ABC")
	sample_3_folder_breseq = checkdir(sample_3_folder / "breseq_output")
	sample_3_folder_output = checkdir(sample_3_folder_breseq / "output")
	sample_3_folder_data = checkdir(sample_3_folder_breseq / "data")
	
	return parent_folder


