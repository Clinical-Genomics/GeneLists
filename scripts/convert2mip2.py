#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function
import sys

def main(argv):
    f = open(argv[0])

    lines = (line.strip() for line in f.readlines())

    mip2_header = ['#Chromosome', 'Gene_start', 'Gene_stop', 'HGNC_symbol', 'Ensembl_gene_id', 'Reduced_penetrance', 'Clinical_db_gene_annotation', 'Disease_associated_transcript']
    mip1_header = ['#Chromosome', 'Gene_start', 'Gene_stop', 'HGNC_symbol', 'Disease_group_pathway', 'Protein_name', 'Symptoms', 'Biochemistry', 'Imaging', 'Disease_trivial_name', 'Trivial_name_short', 'Genetic_disease_model', 'OMIM_gene', 'OMIM_morbid', 'Gene_locus', 'Clinical_db_genome_build', 'UniProt_id', 'Ensembl_gene_id', 'Ensemble_transcript_ID', 'Reduced_penetrance', 'Clinical_db_gene_annotation', 'Disease_associated_transcript']

    print('\t'.join(mip2_header))

    for line in lines:
        # skip comments
        if line.startswith('#'):
            continue

        record = dict(zip(mip1_header, line.split('\t')))
        print('\t'.join([ record[head] for head in mip2_header ]))

if __name__ == '__main__':
    main(sys.argv[1:])
