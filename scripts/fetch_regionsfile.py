#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function
import os
import sys
import re
import argparse
from time import sleep, strftime

from .omim import OMIM
from .ensembl import Ensembl
from .git import getgittag

header = ['Chromosome', 'Gene_start', 'Gene_stop', 'Ensembl_gene_id', 'HGNC_symbol', 'Phenotypic_disease_model', 'OMIM_morbid', 'Ensembl_transcript_to_refseq_transcript', 'Gene_description']

def fill_line(row):
    """Formats a line in the regions file with information found in the row

    Args:
        row (dict): with following keys: Chromosome, Gene_start, Gene_stop, Ensembl_ID, HGNC_symbol

    Returns: a tab delimited string

    """
    # sanity check
    if row['Gene_start'] > row['Gene_stop']:
        row['Gene_start'], row['Gene_stop'] = row['Gene_stop'], row['Gene_start']

    line = [ str(row.get(column_header, '')) for column_header in header ]
    return '\t'.join(line)

def query_omim(data):
    """Queries OMIM to fill in the inheritance models

    Args:
        data (list of dicts): Inner dict represents a row in a gene list

    Yields:
        dict: with the added HGNC symbol prepended to the HGNC_symbol column.
    """

    omim = OMIM(api_key='<fill in key>')
    for line in data:
        if 'HGNC_symbol' in line:
            entry = omim.gene(line['HGNC_symbol'])

            phenotypic_disease_model = omim.parse_phenotypic_disease_model(entry['phenotypes'])

            if phenotypic_disease_model != None:
                line['Phenotypic_disease_model'] = '%s:%s' % (line['HGNC_symbol'], phenotypic_disease_model)

            # add OMIM morbid
            if entry['mim_number'] != None:
                line['OMIM_morbid'] = '%s:%s' % (line['HGNC_symbol'], entry['mim_number'])

            sleep(0.25) # wait for 250ms as according to OMIM specs
        yield line

def get_lines(transcripts, repodir):
    """Generate the regions file, line by line, including header

    Args:
        transcripts (list): list of dicts. Dict keys are equal to the header keys.

    Yield: all lines, including header, of the regions file
    """

    version = getgittag(repodir)
    mod_date = strftime('%Y%m%d')

    # yield some headers
    yield '##Database=<ID=cust000-Research.txt,Version=%s,Date=%s,Acronym=Research,Complete_name=Research,Clinical_db_genome_build=GRCh37.p13' % (version, mod_date)
    yield '#' + '\t'.join(header)
    for transcript in transcripts:
        yield fill_line(transcript)

def main(argv):
    parser = argparse.ArgumentParser(description='Queries EnsEMBL to retrieve all protein coding genes, their transcripts and RefSeqIDs. Outputs a research gene list.')
    parser.add_argument('repodir', default=None, help='The path to the git repo where the research list will be stored. Used to retrieve tag/version number')
    args = parser.parse_args(argv)

    with Ensembl() as ensembldb:
        transcripts = \
            query_omim(
            ensembldb.query_transcripts()
        )

    for line in get_lines(transcripts, args.repodir):
        print(line)

if __name__ == '__main__':
    main(sys.argv[1:])
