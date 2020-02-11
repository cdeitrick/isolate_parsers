from pathlib import Path
from typing import *
from tqdm import tqdm

def checkdir(path: Path) -> Path:
	path = Path(path)
	if not path.exists():
		path.mkdir()
	return path


def define_filenames(folder: Path) -> Dict[str, Path]:
	""" Locates all relevant breseq files."""

	files = {
		'index':       folder / "output" / "index.html",
		'gd':          folder / "output" / "evidence" / "annotated.gd",
		'vcf':         folder / "data" / "output.vcf",
		'summary':     folder / "data" / "summary.json",
		'summaryHtml': folder / "output" / "summary.html",
		'marginal':    folder / "output" / "marginal.html"
	}

	return files


def make_breseq_folders(parent: Path) -> Path:
	checkdir(parent)
	checkdir(parent / "data")
	checkdir(parent / "output")
	checkdir(parent / "output" / "evidence")

	return parent


def move_breseq_folder(source_folder: Path, destination_folder: Path):
	"""
		Moves the relevant breseq files into a new folder. Does not move any of the unneccessary files.
	Parameters
	----------
	source_folder: Source breseq folder.
	destination_folder: The parent folder where the files should be saved. The contents of the source folder
		will be added to a subfolder named after the linked sample.
	"""
	make_breseq_folders(destination_folder)
	files_source = define_filenames(source_folder)
	files_destination = define_filenames(destination_folder)

	for key in files_source.keys():
		filename_source = files_source[key]
		if filename_source.exists():
			filename_destination = files_destination[key]
			filename_destination.write_bytes(filename_source.read_bytes())


def collect_breseq_folders(folder: Path) -> Dict[str, Path]:
	folders = dict()

	for subfolder in folder.iterdir():
		if subfolder.is_file(): continue
		sample_id = subfolder.name
		folders[sample_id] = subfolder

	return folders


def get_parser(args = None):
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

	if args:
		args = parser.parse_args(args)
	else:
		args = parser.parse_args()
	return args


def main(args):
	destination_folder = checkdir(args.destination)
	breseq_folders = collect_breseq_folders(args.path)

	for sample_id, source_folder in tqdm(list(breseq_folders.items())):
		destination_folder_sample = checkdir(destination_folder/ sample_id)
		move_breseq_folder(source_folder, destination_folder_sample)


if __name__ == "__main__":
	import argparse

	debug_args = [
		"/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/isolatparserdata/cefepime",
		"/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/isolatparserdata/cefepime2"
	]
	args = get_parser(debug_args)

	main(args)
