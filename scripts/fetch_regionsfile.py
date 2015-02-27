#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function
import os
import sys
import re
import argparse
from time import sleep, strftime

from .ensembl import Ensembl
from .git import getgittag

header = ['Chromosome', 'Gene_start', 'Gene_stop', 'Ensembl_gene_id', 'HGNC_symbol', 'Ensembl_transcript_to_refseq_transcript', 'Gene_description']

def fill_line(row):
    """Formats a line in the regions file with information found in the row

    Args:
        row (dict): with following keys: Chromosome, Gene_start, Gene_stop, Ensembl_ID, HGNC_symbol

    Returns: a tab delimited string

    """
    # sanity check
    if row['Gene_start'] > row['Gene_stop']:
        row['Gene_start'], row['Gene_stop'] = row['Gene_stop'], row['Gene_start']

    # set a default value for None
    for key,value in row.items():
        row[key] = '' if value == None else str(value)

    line = [ row[key] for key in header ]
    return '\t'.join(line)

def fill(data):
    """Removes #NA's and None's

    Args:
        data (list of dicts): representing the lines and columns in a gene list
    Yields:
        dict: with all missing columns filled in with ''

    """
    defaults=dict((column_name, '') for column_name in header)
    for line in data:
        d = defaults.copy()
        d.update(line)
        yield d

def get_lines(transcripts, repodir):
    """Generate the regions file, line by line, including header

    Args:
        transcripts (list): list of dicts. Dict keys are equal to the header keys.

    Yield: all lines, including header, of the regions file
    """

    version = getgittag(repodir)
    mod_date = strftime('%Y%m%d')

    # yield some headers
    yield '##Database=<ID=cust000-Research.txt,Version=%s,Date=%s,Acronym=Research,Clinical_db_genome_build=GRCh37.p13' % (version, mod_date)
    yield '#' + '\t'.join(header)
    for transcript in transcripts:
        yield fill_line(transcript)

def main(argv):
    parser = argparse.ArgumentParser(description='Queries EnsEMBL to retrieve all protein coding genes, their transcripts and RefSeqIDs. Outputs a research gene list.')
    parser.add_argument('repodir', default=None, help='The path to the git repo where the research list will be stored. Used to retrieve tag/version number')
    args = parser.parse_args(argv)

    with Ensembl() as ensembldb:
        transcripts = \
            fill(
            ensembldb.query_transcripts()
        )

    for line in get_lines(transcripts, args.repodir):
        print(line)

if __name__ == '__main__':
    main(sys.argv[1:])
