from setuptools import setup
from pathlib import Path

Version = "0.1.2"
FOLDER = Path(__file__).parent
README = FOLDER / "README.md"

with README.open() as readmefile:
	LONG_DESCRIPTION = readmefile.read()

setup(
	name = 'isolateparser',
	version = Version,
	packages = [
		'isolateparser',
		'isolateparser.generate',
		'isolateparser.resultparser',
		'isolateparser.resultparser.parsers'
	],

	extras_require = {
		'Additional support for parsing files': ['beautifulsoup4'],
		'To run tests':                         ["pytest"]
	},
	provides = 'isolateparser',
	url = 'https://github.com/cdeitrick/isolate_parsers',
	license = 'MIT',
	author = 'chris deitrick',
	author_email = 'chrisdeitrick1@gmail.com',
	description = 'A set of scripts to convert multiple breseq analyses together and highlight variabls of interest.',
	long_description = LONG_DESCRIPTION,
	long_description_content_type = 'text/markdown',
	install_requires = [
		'pandas>=0.24.0', 'loguru', 'xlrd', 'unidecode'
	],
	tests_requires = ['pytest'],
	classifiers = [
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
	],
	scripts = ["breseqparser"]
)
