from pathlib import Path
from typing import List, Dict, Tuple
import pandas
import math
from breseqparser.file_parsers.parse_index import parse_index_file


def convert_to_muller_format(filenames: List[Path], timepoints: List[int], population_name = "population"):
	tables = list()
	for filename, timepoint in zip(filenames, timepoints):
		variant_table, *_ = parse_index_file(population_name, filename, set_index = False)
		variant_table['timepoint'] = timepoint
		variant_table['freq'] = [float(i[:-1]) for i in variant_table['freq'].values]
		tables.append(variant_table)

	df: pandas.DataFrame = pandas.concat(tables)
	df = df.set_index(['seq id', 'position', 'timepoint'])
	muller_table = dict()
	for index, row in df.iterrows():
		seq, position, timepoint = index
		frequency = row['freq']
		if (seq, position) not in muller_table:
			data = {
				'seq id': seq,
				'position': position,
				timepoint: frequency
			}
			muller_table[seq, position] = data
		else:
			muller_table[seq, position][timepoint] = frequency
	muller_table = pandas.DataFrame(list(muller_table.values()))
	muller_table = muller_table.fillna(0)
	muller_table = muller_table.sort_values(by = ['seq id', 'position'])
	print(muller_table.to_string())

if __name__ == "__main__":
	filenames = [
		Path("/home/cld100/Documents/github/isolate_parsers/tests/data/Population_Output/2K/index.html"),
		Path("/home/cld100/Documents/github/isolate_parsers/tests/data/Population_Output/5K/index.html"),
		Path("/home/cld100/Documents/github/isolate_parsers/tests/data/Population_Output/10K/index.html"),
		Path("/home/cld100/Documents/github/isolate_parsers/tests/data/Population_Output/15K/index.html"),
		Path("/home/cld100/Documents/github/isolate_parsers/tests/data/Population_Output/20K/index.html")
	]
	timepoints = [2, 5, 10, 15, 20]

	convert_to_muller_format(filenames, timepoints)
