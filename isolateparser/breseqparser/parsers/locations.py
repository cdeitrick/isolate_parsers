""" Implements the logic to find specific files in a breseq folder."""
from pathlib import Path
from typing import *

from loguru import logger

FILENAME_INDEX = "index.html"
FILENAME_GD = "annotated.vcf"
FILENAME_VCF = "output.vcf"
FILENAME_SUMMARY = "summary.json"


def get_candidate_paths(parent: Path, name: str) -> List[Path]:
	""" Combines all of the candidate paths for each file for simplicity. The previous method of
		defining candidates individually was getting complicated.
	"""

	candidates = [
		# The usual locations for the index file.
		parent / "breseq" / "output" / name,
		parent / "breseq" / name,
		parent / "breseq_output" / name,
		parent / "breseq_output" / "output" / name,
		parent / "output" / name,
		# The usual locations for the gd file.
		parent / "output" / "evidence" / name,
		parent / "evidence" / name,
		# The VCF and summary.json files
		parent / "breseq_output" / "data" / name,
		parent / "breseq" / "data" / name,
		parent / "data" / name,
		# Other
		parent / name,
		Path(name)
	]
	return candidates


def get_filename(parent: Path, name: str) -> Optional[Path]:
	candidates = get_candidate_paths(parent, name)
	candidate = filesearch(candidates, parent, name)
	if candidate:
		candidate = candidate.absolute()
	return candidate


def get_file_locations(folder: Path) -> Dict[str, Path]:
	result = {
		'index':   get_filename(folder, 'index.html'),
		'gd':      get_filename(folder, 'annotated.gd'),
		'vcf':     get_filename(folder, 'output.vcf'),
		'summary': get_filename(folder, 'summary.json')
	}
	return result


def is_folder_breseq(parent: Path) -> bool:
	"""
		Tests whether the folder is a breseq run folder.
	"""
	indexpath = get_filename(parent, FILENAME_INDEX)
	return indexpath and indexpath.exists()


def filesearch(candidates: Iterable[Path], parent: Path = None, name: str = None) -> Optional[Path]:
	""" Returns the first filename in `candidates` that exists.
		Parameters
		----------
		candidates: List[Path]
			The filenames to search through.
		parent: Path
			Should only be passed alongside `name`
		name: str; defualt None
			If given, will perform a globular search and attempt to find a filename with the name `name`. Include the extension.
	"""
	if (parent and name) and parent.name == name:
		return parent
	candidate = searchcandidates(candidates)
	if candidate is None and parent is not None and name is not None:
		candidates = parent.glob(f"**/{name}")
		candidate = filesearch(candidates)
	return candidate


def searchcandidates(candidates: Iterable[Path]) -> Optional[Path]:
	""" Searches for the first filename in `candidates` which exists. Returns `None` if none of the paths are valid.
		Parameters
		----------
		candidates: List[Path]
			The filenames to test.
	"""
	for candidate in candidates:
		if candidate.exists():
			return candidate
