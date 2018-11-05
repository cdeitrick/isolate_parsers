from pathlib import Path
import pandas

def generate_reference_sequence(snp_table: pandas.DataFrame)->pandas.Series:
	"""
		Generates a fasta file from the full snp table for a breseq run.
	Parameters
	----------
	snp_table: pandas.DataFrame
	- `sampleName`
	- `seq id`
	- `position`
	- `ref`
	- `alt`
	"""
	snp_table = snp_table.reset_index()
	snp_table = snp_table
	mutation_keys = zip(snp_table['seq id'].tolist(), snp_table['position'].tolist())
	unique_keys = sorted(set(mutation_keys))


	reference = snp_table.set_index(['seq id', 'position'])['ref']
	reflist = list()
	for unique_key in unique_keys:
		a, b = unique_key
		value = reference[unique_key]
		if isinstance(value, pandas.DataFrame):
			value = value.iloc[0]
		if isinstance(value, pandas.Series):
			value = value.tolist().pop()
		if isinstance(value, float):
			value = '.'
		refrow = {
			'seq id': a,
			'position': b,
			'ref': value
		}
		if not isinstance(refrow['ref'], str):
			print(unique_key, type(refrow['ref']))
			print(refrow['ref'])
		reflist.append(refrow)

	df = pandas.DataFrame(reflist)
	df = df.set_index(['seq id', 'position'])
	se = df['ref']
	se.name = 'reference'
	return se


def generate_fasta(snp_table):
	snp_table = snp_table.reset_index()
	snp_table = snp_table[snp_table['variantType'] == 'snp']
	reference_sequence:pandas.Series = generate_reference_sequence(snp_table)

	sample_list = snp_table['sampleName'].unique()
	sample_alts = [reference_sequence]
	groups = snp_table.groupby(by = 'sampleName')
	for sample_name, group in groups:
		group = group.set_index(['seq id', 'position'])
		sample_alt:pandas.Series = group['alt']
		sample_alt.name = sample_name

		sample_alt,_ = sample_alt.align(reference_sequence, join = 'right')
		sample_alt = sample_alt.where(sample_alt.notna(), other = reference_sequence)
		sample_alts.append(sample_alt)

	df: pandas.DataFrame = pandas.concat(sample_alts, axis = 1)
	df = df.transpose()
	return df


if __name__ == "__main__":
	from breseqset_parser import parse_breseqset
	_folder = Path("/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/lipuma/pipeline_output/")
	whitelist = Path("/home/cld100/Documents/projects/lipuma/siblingpairA.txt")
	whitelist = whitelist.read_text().split('\n')
	snp_table = parse_breseqset(_folder)

	snp_table = snp_table[snp_table['Sample'].isin(whitelist)]
	#snp_table = snp_table[snp_table['quality'] >= 30]
	df = generate_fasta(snp_table)

	with Path(__file__).with_name("df.unfiltered.fasta").open('w') as file1:
		for index, row in df.iterrows():
			print(index, len(row))
			seq = "".join(row.tolist())
			file1.write(f">{index}\n{seq}\n")

