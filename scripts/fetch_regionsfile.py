#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function
import os
import sys
import re

from .ensembl import Ensembl

def fill_line(row):
    """Formats a line in the regions file with information found in the row

    Args:
        row (dict): with following keys: Chromosome, Gene_start, Gene_stop, Ensembl_ID, HGNC_symbol

    Returns: a tab delimited string

    """

    # sanity check
    Gene_start = row['Gene_start']
    Gene_stop  = row['Gene_stop']
    if row['Gene_start'] > row['Gene_stop']:
        Gene_start, Gene_stop = Gene_stop, Gene_start

    # set a default value for None
    for key,value in row.items():
        row[key] = '' if value == None else str(value)

    return '\t'.join([row['Chromosome'], row['Gene_start'], row['Gene_stop'], row['Ensembl_gene_id'], row['HGNC_symbol'], row['Ensembl_transcript_to_refseq_transcript'], row['Gene_description']])

def main(argv):

    with Ensembl() as ensembldb:
        transcripts = ensembldb.query_transcripts()

    print('#Chromosome	Gene_start	Gene_stop	Ensembl_gene_id	HGNC_symbol	Ensembl_transcript_to_refseq_transcript	Gene_description')
    for transcript in transcripts:
        print(fill_line(transcript))

if __name__ == '__main__':
    main(sys.argv[1:])
