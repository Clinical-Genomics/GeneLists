#!/usr/bin/env python
# encoding: utf-8

import unittest
from time import strftime
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

    def test_query_omim(self):
        """Tests processed results from OMIM
        """
        input_data = [
            { 'HGNC_symbol': 'TPTE' },
            { 'HGNC_symbol': 'ALG9' },
            { 'HGNC_symbol': 'AARS2' }
        ]

        expected_output_lines = [
            {'HGNC_symbol': 'TPTE', 'OMIM_morbid': 'TPTE:604336'},
            {'Phenotypic_disease_model': 'ALG9:608776', 'HGNC_symbol': 'ALG9', 'OMIM_morbid': 'ALG9:606941'},
            {'Phenotypic_disease_model': 'AARS2:614096>AR|615889>AR', 'HGNC_symbol': 'AARS2', 'OMIM_morbid': 'AARS2:612035'}
        ]

        output_lines = [ line for line in fetch_regionsfile.query_omim(input_data) ]
        for i, line in enumerate(output_lines):
            self.assertEqual(line, expected_output_lines[i])

    def test_get_lines(self):
        """Test presence of correctly formatted header and output
        Returns: pass

        """

        version = None
        mod_date = strftime('%Y%m%d')

        expected_output_lines = [
            '##Database=<ID=cust000-Research.txt,Version=%s,Date=%s,Acronym=Research,Clinical_db_genome_build=GRCh37.p13' % (version, mod_date),
            '#Chromosome	Gene_start	Gene_stop	Ensembl_gene_id	HGNC_symbol	Phenotypic_disease_model	OMIM_morbid	Ensembl_transcript_to_refseq_transcript	Gene_description',
            '21	10906201	11029719	ENSG00000166157	TPTE		TPTE:604336	ENSG00000166157:ENST00000298232>NM_199259|ENST00000328758|ENST00000342420>NM_199260|ENST00000361285>NM_199261|ENST00000415664|ENST00000447568	transmembrane_phosphatase_with_tensin_homology',
            '11	111652919	111742305	ENSG00000086848	ALG9	ALG9:608776	ALG9:606941	ENSG00000086848:ENST00000398006>NM_001077690/NM_001077691/NM_001077692/XM_005277725|ENST00000524386|ENST00000524457|ENST00000524671|ENST00000525910|ENST00000526272|ENST00000527212|ENST00000527228|ENST00000527294|ENST00000527714|ENST00000527883|ENST00000529754|ENST00000530723|ENST00000530851|ENST00000531154>NM_024740|ENST00000532374|ENST00000532425	ALG9__alpha-1_2-mannosyltransferase'
        ]

        with Ensembl() as ensembldb:
            transcripts = fetch_regionsfile.query_omim([ ensembldb.query_transcripts(ensembl_id) for ensembl_id in self.data])
            lines = [ line for line in fetch_regionsfile.get_lines(transcripts, '/') ]
            self.assertEqual(lines, expected_output_lines)

if __name__ == '__main__':
    unittest.main()
