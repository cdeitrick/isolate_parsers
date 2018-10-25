from pathlib import Path
from typing import Optional
import pandas
import vcf


def _get_vcf_path(path: Path) -> Path:
	if path.is_dir():
		result = path / "data" / "output.vcf"
		if not result.exists():
			candidates = list(path.glob("**/output.vcf"))
			if len(candidates) != 1:
				message = f"Invalid vcf file path: {path}"
				raise FileNotFoundError(message)
			result = candidates[0]
	elif path.suffix == '.vcf':
		result = path
	else:
		message = f"Invalid VCF path: {path}"
		raise ValueError(message)
	return result


def parse_vcf(path: Path, name: Optional[str] = None):
	"""
		Converts the VCF file generated by breseq into a pandas Dataframe.
	Parameters
	----------
	path:Path
		Either a folder containing a single breseq run or a path to the vcf file itself.
		name: Optional[str]
			The name of the sample will be added as the `sampleName` column if provided.
	"""
	filename = _get_vcf_path(path)
	table = list()
	with filename.open('r') as file1:
		vcf_reader = vcf.Reader(file1)

		for record in vcf_reader:
			alt = record.ALT
			ref = record.REF
			qual = record.QUAL
			depth = record.INFO.get('DP')
			row = {
				'chromosome': record.CHROM,
				'position':   record.POS,
				'alt':        alt,
				'ref':        ref,
				'quality':    qual,
				'readDepth':  depth
			}
			table.append(row)
	df = pandas.DataFrame(table)
	if name:
		df['sampleName'] = name
	return df


if __name__ == "__main__":
	data_folder = Path(__file__).with_name('data')

	d = parse_vcf(data_folder)
	print(d)
