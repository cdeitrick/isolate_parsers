from pathlib import Path
from bs4 import BeautifulSoup
import json
from pprint import pprint
if __name__ == "__main__":
	filename = Path("/home/cld100/Documents/github/isolate_parsers/tests/data/summary_files/summary.json")

	data = json.loads(filename.read_text())

	reads = data['reads']

	total_bases = reads['total_bases']
	total_reads = reads['total_reads']

	total_aligned_bases = reads['total_aligned_bases']
	total_aligned_reads = reads['total_aligned_reads']

	fraction_aligned_bases = reads['total_fraction_aligned_bases']
	fraction_aligned_reads = reads['total_fraction_aligned_reads']

	row = {
		'totalBases': total_bases,
		'totalReads': total_reads,
		'alignedBases': total_aligned_bases,
		'alignedReads': total_aligned_reads,
		'fractionAlignedBases': fraction_aligned_bases,
		'fractionAlignedReads': fraction_aligned_reads
	}


