#!/usr/bin/env python
# encoding: utf-8

import unittest
import pymysql
from scripts import fetch_info

# extra tests
# test removal autosomal models X-linked gene: NDUFA1  300078  Xq23.3  X
# test removal autosomal models Mitochondrial: SURF1   185620  9q23.3  9 
# test removal unchecked and susceptible models: TACO1   612958  17q23.3 17 

class TestFetchInfo(unittest.TestCase):

    def setUp(self):
        fetch_info.conn = pymysql.connect(host='ensembldb.ensembl.org', port=5306, user='anonymous', db='homo_sapiens_core_75_37')

    def test_query_transcripts(self):
        input_data = [
                {'Chromosome': '21', 'Gene_start': 10906201, 'HGNC_symbol': 'TPTE', 'Ensembl_gene_id': 'ENSG00000166157', 'Gene_stop': 11029719, '#Chromosome': '21'},
                {'Chromosome': '11', 'HGNC_symbol': 'ALG9', 'Ensembl_gene_id': 'ENSG00000086848', '#Chromosome': '11', 'Gene_start': 111652919, 'Gene_stop': 111742305}
        ]

        expected_output_lines = [
                {'#Chromosome': '21', 'Ensembl_gene_id': 'ENSG00000166157', 'Ensembl_transcript_to_refseq_transcript': 'ENSG00000166157:ENST00000298232>NM_199259|ENST00000328758|ENST00000342420>NM_199260|ENST00000361285>NM_199261|ENST00000415664|ENST00000447568', 'Chromosome': '21', 'Gene_description': 'transmembrane_phosphatase_with_tensin_homology', 'HGNC_symbol': 'TPTE', 'Gene_start': 10906201, 'Gene_stop': 11029719},
                {'HGNC_symbol': 'ALG9', 'Ensembl_gene_id': 'ENSG00000086848', 'Gene_description': 'ALG9__alpha-1_2-mannosyltransferase', 'Gene_start': 111652919, '#Chromosome': '11', 'Ensembl_transcript_to_refseq_transcript': 'ENSG00000086848:ENST00000398006>NM_001077690/NM_001077691/NM_001077692/XM_005277725|ENST00000524386|ENST00000524457|ENST00000524671|ENST00000525910|ENST00000526272|ENST00000527212|ENST00000527228|ENST00000527294|ENST00000527714|ENST00000527883|ENST00000529754|ENST00000530723|ENST00000530851|ENST00000531154>NM_024740|ENST00000532374|ENST00000532425', 'Chromosome': '11', 'Gene_stop': 111742305}
        ]

        output_lines = [ line for line in fetch_info.query_transcripts(input_data) ]
        for i, line in enumerate(output_lines):
            self.assertEqual(line, expected_output_lines[i])

    def test_query_omim(self):
        """Tests processed results from OMIM
        """
        input_data = [
            { 'HGNC_symbol': 'TPTE' },
            { 'HGNC_symbol': 'ALG9' },
            { 'HGNC_symbol': 'AARS2' }
        ]

        expected_output_lines = [
            {'HGNC_symbol': 'TPTE', 'OMIM_morbid': 'TPTE:604336', 'Gene_locus': '21p11.1'},
            {'Phenotypic_disease_model': 'ALG9:608776', 'HGNC_symbol': 'ALG9', 'OMIM_morbid': 'ALG9:606941', 'Gene_locus': '11q23.1'},
            {'Phenotypic_disease_model': 'AARS2:614096>AR|615889>AR', 'HGNC_symbol': 'AARS2', 'OMIM_morbid': 'AARS2:612035', 'Gene_locus': '6p21.1'}
        ]

        output_lines = [ line for line in fetch_info.query_omim(input_data) ]
        for i, line in enumerate(output_lines):
            self.assertEqual(line, expected_output_lines[i])

if __name__ == '__main__':
    unittest.main()
