#!/usr/bin/env python
# encoding: utf-8

import unittest
from scripts import fetch_regionsfile
from scripts.ensembl import Ensembl

class TestFetchRegionsfile(unittest.TestCase):

    def setUp(self):
        self.data = [ 'ENSG00000166157', 'ENSG00000086848' ]

    def test_fill_line(self):
        expected_output_lines = [
                {'Ensembl_gene_id': 'ENSG00000166157', 'Ensembl_transcript_to_refseq_transcript': 'ENSG00000166157:ENST00000298232>NM_199259|ENST00000328758|ENST00000342420>NM_199260|ENST00000361285>NM_199261|ENST00000415664|ENST00000447568', 'Chromosome': '21', 'Gene_description': 'transmembrane_phosphatase_with_tensin_homology', 'HGNC_symbol': 'TPTE', 'Gene_start': 10906201, 'Gene_stop': 11029719},
                {'HGNC_symbol': 'ALG9', 'Ensembl_gene_id': 'ENSG00000086848', 'Gene_description': 'ALG9__alpha-1_2-mannosyltransferase', 'Gene_start': 111652919, 'Ensembl_transcript_to_refseq_transcript': 'ENSG00000086848:ENST00000398006>NM_001077690/NM_001077691/NM_001077692/XM_005277725|ENST00000524386|ENST00000524457|ENST00000524671|ENST00000525910|ENST00000526272|ENST00000527212|ENST00000527228|ENST00000527294|ENST00000527714|ENST00000527883|ENST00000529754|ENST00000530723|ENST00000530851|ENST00000531154>NM_024740|ENST00000532374|ENST00000532425', 'Chromosome': '11', 'Gene_stop': 111742305}

        ]

        with Ensembl() as ensembldb:
            for i, ensembl_id in enumerate(self.data):
                output_data = ensembldb.query_transcripts(ensembl_id)
                self.assertEqual(output_data, expected_output_lines[i])

if __name__ == '__main__':
    unittest.main()
