from . import locations
from .parse_index import IndexParser
from .parse_gd import GDParser
from .parse_vcf import parse_vcf_file
from .parse_summary import parse_summary_file
from .gdtotable import GDToTable

field_groups = {
	'primary': ['category_id', 'evidence_id', 'parent_ids'],
	'coverage': [],
	'omit': ['key', 'html_gene_name'],
	'codon': ['codon_new_seq', 'codon_number', 'codon_position', ''],
	'amino': ['aa_new_seq', 'aa_position', 'aa_ref_seq'],
	'gene': ['locus_tag', 'locus_tags_overlapping', 'gene_name', 'gene_position', 'gene_product', 'gene_strand', 'genes_overlapping'],
	'mutation': ['mutation_category', 'new_seq', 'seq_id', 'position', 'snp_type']
}

snp_columns = {
    "aa_new_seq": "V",
    "aa_position": "63",
    "aa_ref_seq": "A",
    "category_id": "SNP",
    "codon_new_seq": "GTC",
    "codon_number": "63",
    "codon_position": "2",
    "codon_ref_seq": "GCC",
    "evidence_id": "79",
    "gene_name": "JNAKOBFD_06328",
    "gene_position": "188",
    "gene_product": "[locus_tag=Bcen2424_3247] [db_xref=InterPro:IPR000089,InterPro:IPR000103, InterPro:IPR000205,InterPro:IPR000815,InterPro:IPR001100, InterPro:IPR001327,InterPro:IPR002218,InterPro:IPR003016, InterPro:IPR004099,InterPro:IPR006258,InterPro:IPR012999, InterPro:IPR013027][protein=dihydrolipoamide dehydrogenase] [protein_id=ABK09992.1][location=86394..88181] [gbkey=CDS]",
    "gene_strand": ">",
    "genes_overlapping": "JNAKOBFD_06328",
    "html_gene_name": "<i>JNAKOBFD_06328</i>&nbsp;&rarr;",
    "locus_tag": "JNAKOBFD_06328",
    "locus_tags_overlapping": "JNAKOBFD_06328",
    "mutation_category": "snp_nonsynonymous",
    "new_seq": "T",
    "parent_ids": "350",
    "position": "802766",
    "seq_id": "3",
    "snp_type": "nonsynonymous",
    "transl_table": "11"
}

ins_columns = {
    "category_id": "INS",
    "evidence_id": "77",
    "gene_name": "JNAKOBFD_06217",
    "gene_position": "coding (87/273 nt)",
    "gene_product": "hypothetical protein",
    "gene_strand": "<",
    "genes_inactivated": "JNAKOBFD_06217",
    "html_gene_name": "<i>JNAKOBFD_06217</i>&nbsp;&larr;",
    "locus_tag": "JNAKOBFD_06217",
    "locus_tags_inactivated": "JNAKOBFD_06217",
    "mutation_category": "small_indel",
    "new_seq": "G",
    "parent_ids": "340",
    "position": "682313",
    "repeat_length": "1",
    "repeat_new_copies": "7",
    "repeat_ref_copies": "6",
    "repeat_seq": "G",
    "seq_id": "3"
}

del_columns = {
    "category_id": "DEL",
    "evidence_id": "80",
    "gene_name": "JNAKOBFD_06387",
    "gene_product": "JNAKOBFD_06387",
    "genes_inactivated": "JNAKOBFD_06387",
    "html_gene_name": "<i>JNAKOBFD_06387</i>",
    "locus_tag": "[JNAKOBFD_06387]",
    "locus_tags_inactivated": "JNAKOBFD_06387",
    "mutation_category": "large_deletion",
    "parent_ids": "364",
    "position": "1",
    "seq_id": "5",
    "size": "1248"
}
