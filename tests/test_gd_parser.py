import unittest
from pathlib import Path

import pandas
import pytest

from isolateparser.breseqoutputparser import file_parsers


@pytest.fixture
def mutation_snp() -> str:
	MUTATION_SNP = "SNP	3	48	REL606	380188	C	aa_new_seq=L	aa_position=239	aa_ref_seq=F	codon_new_seq=TTG	codon_number=239	" \
				   "codon_position=3	codon_ref_seq=TTT	gene_name=araJ	gene_position=717	gene_product=predicted transporter	gene_strand=<	" \
				   "genes_overlapping=araJ	html_gene_name=<i>araJ</i>&nbsp;&larr;	locus_tag=ECB_00344	locus_tags_overlapping=ECB_00344	" \
				   "mutation_category=snp_nonsynonymous	snp_type=nonsynonymous	transl_table=11"
	return MUTATION_SNP


@pytest.fixture
def mutation_indel() -> str:
	MUTATION_INDEL = "INS	5	51	REL606	475292	G	gene_name=ybaL	gene_position=coding (14/1677 nt)	" \
					 "gene_product=predicted transporter with NAD(P)-binding Rossmann-fold domain	gene_strand=<	genes_inactivated=ybaL	" \
					 "html_gene_name=<i>ybaL</i>&nbsp;&larr;	locus_tag=ECB_00429	locus_tags_inactivated=ECB_00429	mutation_category=small_indel"
	return MUTATION_INDEL


@pytest.fixture
def evidence_ra() -> str:
	EVIDENCE_RA = "RA	81	.	REL606	3037598	0	A	C	aa_new_seq=G	aa_position=466	aa_ref_seq=V	bias_e_value=15085.6	" \
				  "bias_p_value=0.00325836	codon_new_seq=GGT	codon_number=466	codon_position=2	codon_ref_seq=GTT	" \
				  "consensus_reject=FREQUENCY_CUTOFF	consensus_score=47.7	fisher_strand_p_value=0.000365526	frequency=3.243e-01	" \
				  "gene_name=ECB_02838	gene_position=1397	gene_product=GspD, hypothetical type II secretion protein	gene_strand=<	" \
				  "html_gene_name=<i>ECB_02838</i>&nbsp;&larr;	ks_quality_p_value=1	locus_tag=ECB_02838	major_base=A	major_cov=15/10	" \
				  "major_frequency=6.757e-01	minor_base=C	minor_cov=0/12	new_cov=0/12	new_seq=C	polymorphism_frequency=3.243e-01	" \
				  "polymorphism_score=10.9	prediction=polymorphism	ref_cov=15/10	ref_seq=A	snp_type=nonsynonymous	total_cov=15/22	transl_table=11"
	return EVIDENCE_RA


data_folder = Path(__file__).parent / 'data' / 'Clonal_Output' / 'breseq_output'
gd_filename = data_folder / 'output' / 'evidence' / 'annotated.gd'


def test_sort_gd_file_rows():
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
	mutations, evidence = file_parsers.parse_gd._sort_gd_file_rows(string)

	assert truth_mutations == mutations
	assert truth_evidence == evidence


def test_get_row_position_keys_by_mutation():
	row_type = 'snp'
	other = "REL606	161041	G	aa_new_seq=H	aa_position=302	aa_ref_seq=N".split('\t')
	truth = {
		'seqId':    'REL606',
		'position': '161041',
		'new_seq':  'G'
	}
	result = file_parsers.parse_gd._get_row_position_and_sequence(row_type, other)

	assert truth == result


def test_extract_reference_base_from_codon():
	codon = 'TAC'
	position = '1'

	assert 'T' == file_parsers.parse_gd._extract_reference_base_from_codon(codon, position)

	codon = 'GAT'
	position = '3'

	assert 'T' == file_parsers.parse_gd._extract_reference_base_from_codon(codon, position)


def test_parse_annotated_gd_file_row_indel(mutation_indel):
	expected_details = {
		# 'seqId': 'REL606',
		# 'position': '475292',
		'new_seq':                'G',
		'gene_name':              'ybaL',
		'gene_position':          'coding (14/1677 nt)',
		'gene_product':           'predicted transporter with NAD(P)-binding Rossmann-fold domain',
		'gene_strand':            '<',
		'genes_inactivated':      'ybaL',
		'html_gene_name':         '<i>ybaL</i>&nbsp;&larr;',
		'locus_tag':              'ECB_00429',
		'locus_tags_inactivated': 'ECB_00429',
		'mutation_category':      'small_indel'
	}

	expected_mutation = file_parsers.parse_gd.Mutation('ins', '5', '51', seqId = 'REL606', position = 475292, details = expected_details)
	test_mutation = file_parsers.parse_gd.parse_annotated_gd_file_row(mutation_indel.split('\t'))

	assert test_mutation == expected_mutation


def test_parse_annotated_gd_file_row_evidence(evidence_ra):
	expected_details = {
		'seqId':                 'REL606',
		'position':              '3037598',
		'insert_position':       '0',
		'ref_base':              'A',
		'new_base':              'C',
		'aa_new_seq':            'G',
		'aa_position':           '466',
		'aa_ref_seq':            'V',
		'bias_e_value':          '15085.6',
		'bias_p_value':          '0.00325836',
		'codon_new_seq':         'GGT',
		'codon_number':          '466',
		'codon_position':        '2',
		'codon_ref_seq':         'GTT',
		'consensus_reject':      'FREQUENCY_CUTOFF',
		'consensus_score':       '47.7',
		'fisher_strand_p_value': '0.000365526',
		'frequency':             '3.243e-01',
		'gene_name':             'ECB_02838',
		'gene_position':         '1397',
		'gene_product':          'GspD, hypothetical type II secretion protein',
		'gene_strand':           '<',
		'html_gene_name':        '<i>ECB_02838</i>&nbsp;&larr;',
		'ks_quality_p_value':    '1',
		'locus_tag':             'ECB_02838',
		'major_base':            'A',
		'major_cov':             '15/10',
		'major_frequency':       '6.757e-01',
		'minor_base':            'C',
		'minor_cov':             '0/12', 'new_cov': '0/12', 'new_seq': 'C', 'polymorphism_frequency': '3.243e-01', 'polymorphism_score': '10.9',
		'prediction':            'polymorphism', 'ref_cov': '15/10', 'ref_seq': 'A', 'snp_type': 'nonsynonymous', 'total_cov': '15/22',
		'transl_table':          '11'
	}
	expected_evidence = file_parsers.parse_gd.Evidence('ra', '81', '.', details = expected_details)
	truth_evidence = file_parsers.parse_gd.parse_annotated_gd_file_row(evidence_ra.split('\t'))
	assert truth_evidence == expected_evidence


def test_parse_annotated_gd_file_row_snp(mutation_snp):
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
	result = file_parsers.parse_gd.parse_annotated_gd_file_row(mutation_snp.split('\t'))

	assert isinstance(result, file_parsers.parse_gd.Mutation)

	assert 'snp' == result.type
	assert '3' == result.id
	assert '48' == result.parentId
	assert 'REL606' == result.seqId
	assert 380188 == result.position

	assert ('REL606', 380188) == result.key

	assert 'araJ' == result.get('gene_name')

	assert set(details.keys()) == set(result.details.keys())
	assert details == result.details


def test_extract_mutations():
	test_list = [
		file_parsers.parse_gd.Mutation('mtype', 'id', 'pid', 'sequence', 1, {}),
		file_parsers.parse_gd.Evidence('etype', 'id', 'parent', {}),
		file_parsers.parse_gd.Mutation('ttype', 'mid', 'mparent', 'seq', 23, {})
	]
	truth_list = [
		file_parsers.parse_gd.Mutation('mtype', 'id', 'pid', 'sequence', 1, {}),
		file_parsers.parse_gd.Mutation('ttype', 'mid', 'mparent', 'seq', 23, {})
	]

	assert truth_list == file_parsers.parse_gd._extract_mutations(test_list)


def test_parse_annotated_gd_file(mutation_snp, mutation_indel, evidence_ra):
	string = f"{mutation_snp}\n{mutation_indel}\n{evidence_ra}"
	gd_data = file_parsers.parse_gd.parse_annotated_gd_file(string)

	expected_data = [
		file_parsers.parse_gd.parse_annotated_gd_file_row(mutation_snp.split('\t')),
		file_parsers.parse_gd.parse_annotated_gd_file_row(mutation_indel.split('\t')),
		file_parsers.parse_gd.parse_annotated_gd_file_row(evidence_ra.split('\t'))
	]
	assert gd_data == expected_data


def test_evidence_object():
	test_evidence = file_parsers.parse_gd.Evidence(
		'RA', 'id1', 'id2', {'seq_id': 'REL606', 'position': '87812', 'bias_e_value': '8.16284', 'bias_p_value': '1.7631e-06'}
	)
	assert 'RA' == test_evidence.type
	assert '8.16284' == test_evidence.get('bias_e_value')
	assert 'REL606' == test_evidence.seqId
	assert 87812 == test_evidence.position
	assert ('REL606', 87812) == test_evidence.key


def test_generate_mutation_table(mutation_snp, mutation_indel):
	mutations = [
		file_parsers.parse_gd.parse_annotated_gd_file_row(mutation_indel.split('\t')),
		file_parsers.parse_gd.parse_annotated_gd_file_row(mutation_snp.split('\t'))
	]

	expected_columns = [
		'description', 'gene', 'mutation', 'position', 'seq id',
		'baseAlt', 'aminoAlt', 'aminoRef', 'codonAlt', 'codonRef',
		'locusTag', 'mutationCategory'
	]
	expected_table = pandas.DataFrame(
		[file_parsers.parse_gd._convert_mutation_to_dictionary(i) for i in mutations]
	)[expected_columns]
	test_table = file_parsers.parse_gd.generate_mutation_table(mutations)
	pandas.testing.assert_frame_equal(expected_table, test_table)


def test_convert_mutation_to_dictionary(mutation_snp):
	mutation = file_parsers.parse_gd.parse_annotated_gd_file_row(mutation_snp.split('\t'))
	expected_dictionary = {
		'description':      'predicted transporter',
		'gene':             'araJ',
		'mutation':         '',
		'position':         380188,
		'seq id':           'REL606',
		'baseAlt':          'C',
		'aminoRef':         'F',
		'aminoAlt':         'L',
		'codonRef':         'TTT',
		'codonAlt':         'TTG',
		'locusTag':         'ECB_00344',
		'mutationCategory': 'snp_nonsynonymous'
	}
	test_dictionary = file_parsers.parse_gd._convert_mutation_to_dictionary(mutation)
	assert test_dictionary == expected_dictionary


def test_parse_gd_file(mutation_snp, mutation_indel, evidence_ra):
	filename = f"{mutation_snp}\n{mutation_indel}\n{evidence_ra}"
	expected_columns = [
		'description', 'gene', 'mutation', 'position', 'seq id',
		'baseAlt', 'aminoAlt', 'aminoRef', 'codonAlt', 'codonRef',
		'locusTag', 'mutationCategory'
	]

	expected_data = [
		file_parsers.parse_gd.parse_annotated_gd_file_row(mutation_snp.split('\t')),
		file_parsers.parse_gd.parse_annotated_gd_file_row(mutation_indel.split('\t'))
	]
	expected_table = pandas.DataFrame(
		[file_parsers.parse_gd._convert_mutation_to_dictionary(i) for i in expected_data]
	)[expected_columns].set_index(keys = ['seq id', 'position'])

	test_table = file_parsers.parse_gd.parse_gd_file(filename, set_index = True)

	pandas.testing.assert_frame_equal(expected_table, test_table)


if __name__ == "__main__":
	unittest.main()
