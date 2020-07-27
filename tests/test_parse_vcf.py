from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas
import pytest

from isolateparser.resultparser.parsers import parse_vcf

data_folder = Path(__file__).parent / 'data' / 'Clonal_Output' / 'breseq_output'
vcf_filename_real = data_folder / "data" / "output.vcf"


@pytest.fixture
def vcf_filename(tmp_path) -> Path:
	contents = """
		##fileformat=VCFv4.1
		##fileDate
		##source=breseq_GD2VCF_converterter
		##contig=<ID=REL606,length=4629812>
		##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
		##INFO=<ID=AD,Number=1,Type=Float,Description="Allele Depth (avg read count)">
		##INFO=<ID=DP,Number=1,Type=Float,Description="Total Depth (avg read count)">
		##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
		REL606	16971	.	GG	GGGGTGGTTGTACTCA	.	PASS	AF=1.0000;AD=59.5;DP=59.5
		REL606	161041	.	T	G	190.1	PASS	AF=1.0000;AD=61;DP=61
		REL606	380188	.	A	C	98.1	PASS	AF=1.0000;AD=31;DP=32
		REL606	430835	.	C	T	135.0	PASS	AF=1.0000;AD=46;DP=48
		REL606	475292	.	G	GG	94.7	PASS	AF=1.0000;AD=30;DP=32
	"""
	#TODO add the scenario where a deletion is follwoed by a SNP. VCF annotates these as idenitical positions.

	f = tmp_path / "vcf_file.vcf"
	f.write_text("\n".join(i.strip() for i in contents.split('\n')))
	return f


@dataclass
class Record:
	""" Easier to emulate a VCF Record rather than instantiate vcf._Record"""
	ALT: List[str]
	REF: str
	QUAL: Optional[float]
	INFO: Dict[str, Any]
	CHROM: str
	POS: int
	var_type: str


def test_convert_record_to_dictionary():
	record = Record(
		CHROM = 'REL606',
		POS = 16971,
		REF = "GG",
		ALT = [
			"GGGGTGGTTGTACTGACCCCAAAAAGTTG"
		],
		var_type = 'indel',
		INFO = {'AF': [1.0], 'AD': 59.5, 'DP': 59.5},
		QUAL = None
	)
	expected = {
		parse_vcf.VCFColumns.sequence_id:  'REL606',
		parse_vcf.VCFColumns.position:     16971,
		parse_vcf.VCFColumns.alternate:    record.ALT[0],
		parse_vcf.VCFColumns.reference:    "GG",
		parse_vcf.VCFColumns.quality:      None,
		parse_vcf.VCFColumns.depth:        59.5,
		parse_vcf.VCFColumns.variant_type: 'indel'
	}
	output = parse_vcf._convert_record_to_dictionary(record)

	assert output == expected


def test_parse_vcf_file(vcf_filename):
	# `seq id` and `position` are in the index.
	expected_columns = ['alt', 'ref', 'quality', 'readDepth', 'variantType']

	vcf_table = parse_vcf.parse_vcf_file(vcf_filename)
	expected_index = [('REL606', 16971), ("REL606", 161041), ("REL606", 380188), ("REL606", 430835), ("REL606", 475292)]
	assert list(vcf_table.columns) == expected_columns
	assert list(vcf_table.index) == list(expected_index)

	snp_table = vcf_table[vcf_table['variantType'] == 'snp']
	expected_snp_index = [("REL606", 161041), ("REL606", 380188), ("REL606", 430835)]
	assert list(snp_table.index) == expected_snp_index


def test_convert_vcf_to_table(vcf_filename):
	record = Record(
		CHROM = 'REL606',
		POS = 16971,
		REF = "GG",
		ALT = ["GGGGTGGTTGTACTCA"],
		var_type = 'indel',
		INFO = {'AF': [1.0], 'AD': 59.5, 'DP': 59.5},
		QUAL = None
	)
	expected = {
		parse_vcf.VCFColumns.sequence_id:  'REL606',
		parse_vcf.VCFColumns.position:     16971,
		parse_vcf.VCFColumns.alternate:    record.ALT[0],
		parse_vcf.VCFColumns.reference:    "GG",
		parse_vcf.VCFColumns.quality:      None,
		parse_vcf.VCFColumns.depth:        59.5,
		parse_vcf.VCFColumns.variant_type: 'indel'
	}
	table = parse_vcf._convert_vcf_to_table(vcf_filename)

	assert table[0] == expected


