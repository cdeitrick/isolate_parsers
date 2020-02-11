from dataclasses import dataclass
from pathlib import Path
from typing import *

import pytest

from isolateparser.breseqparser.parsers import locations

folder_source = Path(__file__).parent / "data"


@dataclass
class BreseqFolder:
	""" Holds the paths for a Breseq output folder, since Fixtures only return a single object."""
	parent: Path
	samples: List[Dict[str, Path]]


def checkdir(path: Path) -> Path:
	path = Path(path)
	if not path.exists():
		path.mkdir()
	return path


def make_file(folder: Path, name: str) -> Path:
	filename = folder / name
	filename.touch()
	return filename


def make_breseq_folder(folder: Path) -> Dict[str, Path]:
	"""
	Mocks a breseq run for a single sample.
	.folder
	|---- output/
	|----|---- index.html
	|---- data/
	|----|---- output.vcf
	"""
	folder = checkdir(folder)
	folder_output = checkdir(folder / "output")
	folder_data = checkdir(folder / "data")
	folder_evidence = checkdir(folder_output / "evidence")

	filename_index = make_file(folder_output, 'index.html')
	filename_vcf = make_file(folder_data, 'output.vcf')
	filename_gd = make_file(folder_evidence, "annotated.gd")
	filename_summary = make_file(folder_data, "summary.json")

	result = {
		'index':   filename_index,
		'gd':      filename_gd,
		'summary': filename_summary,
		'vcf':     filename_vcf
	}

	return result


@pytest.fixture
def structure1(tmp_path) -> BreseqFolder:
	"""
	.parent
	|---- CeN0L1G5/
	|----|---- data/
	|----|---- output/
	|---- CeN0L1G3/
	|----|---- data/
	|----|---- output/
	|---- CeN0L1G29/
	|----|---- data/
	|----|---- output/
	"""

	parent_folder = checkdir(tmp_path / "cefepime")

	filenames_sample_1 = make_breseq_folder(parent_folder / "CeN0L1G5")
	filenames_sample_2 = make_breseq_folder(parent_folder / "CeN0L1G3")
	filenames_sample_3 = make_breseq_folder(parent_folder / "CeN0L1G29")

	samples = [
		filenames_sample_1,
		filenames_sample_2,
		filenames_sample_3
	]

	result = BreseqFolder(
		parent = parent_folder,
		samples = samples
	)

	return result


@pytest.fixture
def structure2(tmp_path) -> BreseqFolder:
	"""
	.parent
	|---- sample1/
	|----|---- breseq
	|---- sample2/
	|----|---- breseq
	|---- sample3/
	|----|---- breseq
	"""

	parent_folder = checkdir(tmp_path / "cefepime")
	folder_1 = checkdir(parent_folder / "CeN0L1G5")
	folder_2 = checkdir(parent_folder / "CeN0L1G3")
	folder_3 = checkdir(parent_folder / "CeN0L1G29")

	filenames_sample_1 = make_breseq_folder(folder_1 / "breseq")
	filenames_sample_2 = make_breseq_folder(folder_2 / "breseq")
	filenames_sample_3 = make_breseq_folder(folder_3 / "breseq")

	result = BreseqFolder(
		parent = parent_folder,
		samples = [
			filenames_sample_1,
			filenames_sample_2,
			filenames_sample_3
		]
	)

	return result


def test_get_file_locations(structure1):
	sample_1 = locations.get_file_locations(structure1.parent / "CeN0L1G5")
	sample_2 = locations.get_file_locations(structure1.parent / "CeN0L1G3")
	sample_3 = locations.get_file_locations(structure1.parent / "CeN0L1G29")

	expected_sample_1 = structure1.samples[0]
	expected_sample_2 = structure1.samples[1]
	expected_sample_3 = structure1.samples[2]

	assert sample_1 == expected_sample_1
	assert sample_2 == expected_sample_2
	assert sample_3 == expected_sample_3


def test_get_file_locations_index(structure1):
	sample_names = ['CeN0L1G5', 'CeN0L1G3', 'CeN0L1G29']
	for sample_name in sample_names:
		folder = structure1.parent / sample_name
		result = locations.get_filename(folder, "index.html")
		expected = folder / "output" / "index.html"

		assert result == expected


def test_get_file_locations_gd(structure1):
	sample_names = ['CeN0L1G5', 'CeN0L1G3', 'CeN0L1G29']
	for sample_name in sample_names:
		folder = structure1.parent / sample_name
		result = locations.get_filename(folder, 'annotated.gd')
		expected = folder / "output" / "evidence" / "annotated.gd"

		assert result == expected


def test_get_file_locations_vcf(structure1):
	sample_names = ['CeN0L1G5', 'CeN0L1G3', 'CeN0L1G29']
	for sample_name in sample_names:
		folder = structure1.parent / sample_name
		result = locations.get_filename(folder, 'output.vcf')
		expected = folder / "data" / "output.vcf"

		assert result == expected


def test_get_file_locations_summary(structure1):
	sample_names = ['CeN0L1G5', 'CeN0L1G3', 'CeN0L1G29']
	for sample_name in sample_names:
		folder = structure1.parent / sample_name
		result = locations.get_filename(folder, 'summary.json')
		expected = folder / "data" / "summary.json"

		assert result == expected

def test_is_breseq_folder(structure1):
	sample_names = ['CeN0L1G5', 'CeN0L1G3', 'CeN0L1G29']
	for sample_name in sample_names:
		folder = structure1.parent / sample_name
		result = locations.get_filename(folder, 'summary.json')
		assert result