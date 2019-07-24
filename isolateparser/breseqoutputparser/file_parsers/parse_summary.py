from pathlib import Path
import json
from typing import Dict, Optional

def get_summary_filename(folder:Path)->Optional[Path]:
	# Assume the folder is a breseq folder.

	expected = folder / "data" / "sumary.json"
	if not expected.exists():
		try:
			found = list(folder.glob("**/summary.json"))[0]
		except IndexError:
			found = None
	else:
		found = expected
	return found

def parse_summary_file(path:Path, sample_id:str, sample_name:Optional[str] = None)->Dict[str,str]:
	if sample_name is None: sample_name = sample_id
	if path.is_dir():
		filename = path / "data" / "summary.json"
	else:
		filename = path

	data = json.loads(filename.read_text())

	reads = data['reads']

	total_bases = reads['total_bases']
	total_reads = reads['total_reads']

	total_aligned_bases = reads['total_aligned_bases']
	total_aligned_reads = reads['total_aligned_reads']

	fraction_aligned_bases = reads['total_fraction_aligned_bases']
	fraction_aligned_reads = reads['total_fraction_aligned_reads']

	row = {
		'sampleId':				sample_id,
		'sampleName':			sample_name,
		'totalBases':           total_bases,
		'totalReads':           total_reads,
		'alignedBases':         total_aligned_bases,
		'alignedReads':         total_aligned_reads,
		'fractionAlignedBases': fraction_aligned_bases,
		'fractionAlignedReads': fraction_aligned_reads
	}
	return row



