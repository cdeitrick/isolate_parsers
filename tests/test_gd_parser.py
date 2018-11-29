import unittest

from breseqparser.file_parsers.parse_gd import *
# Need to also import private functions
from breseqparser.file_parsers.parse_gd import _extract_reference_base_from_codon, _get_row_position_and_sequence, _sort_gd_file_rows

data_folder = Path(__file__).parent / 'data' / 'Clonal_Output' / 'breseq_output'
gd_filename = data_folder / 'output' / 'evidence' / 'annotated.gd'


class TestGDParser(unittest.TestCase):
	def setUp(self):
		self.maxDiff = 1000

	def test_get_gd_filename(self):
		self.assertEqual(gd_filename, get_gd_filename(gd_filename))

	def test_sort_gd_file_rows(self):
		string = "SNP	12	66	REL606	1286699	A\n" \
				 "SNP	13	68	REL606	1329516	T\n" \
				 "MOB	14	128,129	REL606	1733647	IS150	-1	3\n" \
				 "RA	40	.	REL606	87812	0	T	G	bias_e_value=8.16284\n" \
				 "RA	63	.	REL606	1117181	0\n" \
				 "RA	70	.	REL606	1536774	0	G	T	bias_e_value=1.00619\n"

		truth_mutations = [
			['SNP', '12', '66', 'REL606', '1286699', 'A'],
			['SNP', '13', '68', 'REL606', '1329516', 'T'],
			['MOB', '14', '128,129', 'REL606', '1733647', 'IS150', '-1', '3']
		]
		truth_evidence = [
			['RA', '40', '.', 'REL606', '87812', '0', 'T', 'G', 'bias_e_value=8.16284'],
			['RA', '63', '.', 'REL606', '1117181', '0'],
			['RA', '70', '.', 'REL606', '1536774', '0', 'G', 'T', 'bias_e_value=1.00619']
		]
		mutations, evidence = _sort_gd_file_rows(string)

		self.assertListEqual(truth_mutations, mutations)
		self.assertListEqual(truth_evidence, evidence)

	def test_get_row_position_keys_by_mutation(self):
		row_type = 'snp'
		other = "REL606	161041	G	aa_new_seq=H	aa_position=302	aa_ref_seq=N".split('\t')
		truth = {
			'seqId':    'REL606',
			'position': '161041',
			'new_seq':  'G'
		}
		result = _get_row_position_and_sequence(row_type, other)

		self.assertDictEqual(truth, result)

	def test_extract_reference_base_from_codon(self):
		codon = 'TAC'
		position = '1'

		self.assertEqual('T', _extract_reference_base_from_codon(codon, position))

		codon = 'GAT'
		position = '3'

		self.assertEqual('T', _extract_reference_base_from_codon(codon, position))

	def test_parse_annotated_gd_file_row_snp(self):
		string = "SNP	3	48	REL606	380188	C	aa_new_seq=L	aa_position=239	aa_ref_seq=F	codon_new_seq=TTG	codon_number=239	" \
				 "codon_position=3	codon_ref_seq=TTT	gene_name=araJ	gene_position=717	gene_product=predicted transporter	gene_strand=<	" \
				 "genes_overlapping=araJ	html_gene_name=<i>araJ</i>&nbsp;&larr;	locus_tag=ECB_00344	locus_tags_overlapping=ECB_00344	" \
				 "mutation_category=snp_nonsynonymous	snp_type=nonsynonymous	transl_table=11"
		row = string.split('\t')
		details = {
			'seqId':                  'REL606',
			'position':               '380188',
			'new_seq':                'C',
			'aa_new_seq':             'L',
			'aa_position':            '239',
			'aa_ref_seq':             'F',
			'codon_new_seq':          'TTG',
			'codon_number':           '239',
			'codon_position':         '3',
			'codon_ref_seq':          'TTT',
			'gene_name':              'araJ',
			'gene_position':          '717',
			'gene_product':           'predicted transporter',
			'gene_strand':            '<',
			'genes_overlapping':      'araJ',
			'html_gene_name':         '<i>araJ</i>&nbsp;&larr;',
			'locus_tag':              'ECB_00344',
			'locus_tags_overlapping': 'ECB_00344',
			'mutation_category':      'snp_nonsynonymous',
			'snp_type':               'nonsynonymous',
			'transl_table':           '11'
		}
		result = parse_annotated_gd_file_row(row)

		self.assertIsInstance(result, Mutation)

		self.assertEqual('snp', result.type)
		self.assertEqual('3', result.id)
		self.assertEqual('48', result.parentId)
		self.assertEqual('REL606', result.seqId)
		self.assertEqual(380188, result.position)

		self.assertEqual(('REL606', 380188), result.key)

		self.assertEqual('araJ', result.get('gene_name'))
		self.assertIsNone(result.get('asdasdsafdas'))

		self.assertSetEqual(set(details.keys()), set(result.details.keys()))
		self.assertDictEqual(details, result.details)


if __name__ == "__main__":
	unittest.main()
