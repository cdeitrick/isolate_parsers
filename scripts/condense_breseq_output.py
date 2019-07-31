import shutil
from pathlib import Path
from typing import List, Optional, Tuple

from loguru import logger


def get_file(folder: Path, name: str) -> Optional[Path]:
	""" Attempts to retrieve the specified file using `name`"""
	candidates = list(folder.glob(f"**/{name}"))

	try:
		return candidates[0]
	except IndexError:
		return None


def select_files(folder: Path) -> List[Path]:
	""" Locates all relevant breseq files."""
	filename_index = get_file(folder, "index.html")
	filename_marginal = get_file(folder, "marginal.html")
	filename_summary_html = get_file(folder, "summary.html")
	filename_summary_json = get_file(folder, "summary.json")
	filename_gd = get_file(folder, "annotated.gd")
	filename_vcf = get_file(folder, "output.vcf")

	files = [
		filename_index, filename_marginal, filename_summary_html, filename_summary_json, filename_gd, filename_vcf
	]

	return files


def move_breseq(sample_id: str, source_folder: Path, destination_folder: Path):
	"""
		Moves the relevant breseq files into a new folder. Does not move any of the unneccessary files.
	Parameters
	----------
	sample_id:str
	source_folder: Source breseq folder.
	destination_folder: The parent folder where the files should be saved. The contents of the source folder
		will be added to a subfolder named after the linked sample.
	"""
	breseq_files = select_files(source_folder)

	sample_destination = destination_folder / sample_id
	if not sample_destination.exists():
		sample_destination.mkdir()

	for filename in breseq_files:
		destination:Path = sample_destination / filename.name
		shutil.copyfile(filename, destination)


def move_breseq_folder(folders: List[Tuple[str, Path]], destination_folder: Path):
	""" Move a list of breseq folders. The list should be a tuple of sampleIds and breseq folder locations."""
	for sample_id, subfolder in folders:
		move_breseq(sample_id, subfolder, destination_folder)


def collect_breseq_folders(folder: Path) -> List[Tuple[str, Path]]:
	folders = list()

	for subfolder in folder.iterdir():
		if subfolder.is_file(): continue
		sample_id = subfolder.name
		breseq_folder = list(subfolder.glob("**/breseq/"))
		try:
			breseq_folder = breseq_folder[0]
			folders.append((sample_id, breseq_folder))
		except IndexError:
			logger.warning(f"Cannot locate a breseq folder in {subfolder}")
			continue

	return folders


if __name__ == "__main__":
	import argparse

	parser = argparse.ArgumentParser()

	parser.add_argument(
		"path",
		help = "Either a breseq folder or a folder of breseq folders.",
		type = Path
	)

	parser.add_argument(
		"destination",
		help = "THe parent folder to save the breseq files to.",
		type = Path
	)

	parser.add_argument(
		"--sample_id",
		help = "The sample Id, if only moving a single folder",
		type = str
	)
	debug_args = [
		"/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/lipuma/pipelines/SC1360/",
		"/home/cld100/Documents/sandbox/"
	]
	args = parser.parse_args(debug_args)

	breseq_folders = collect_breseq_folders(args.path)

	move_breseq_folder(breseq_folders, args.destination)
