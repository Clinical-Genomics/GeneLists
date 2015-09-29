#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function
import sys

from .omim import OMIM

def list2dict(header, data):
    """Will convert each row in the data from a list to dict using the header list as keys.

    Args:
        header (list): A list containing the keys for the dict generation
        data (list of lists): Inner list represents a row in a gene list
    Yields:
        dict: the next dictified line
    """
    for line in data:
        yield dict(zip(header, line))


def get_raw_omim(data):
    """Retrieves the raw OMIM JSON.

    Args:
        data (list): a list of dicts

    Yields:
        dict: with a raw JSON object include under key 'raw'
    """

    omim = OMIM(api_key='<fill in key>')

    for line in data:
        entry = omim.search_gene_raw(line['HGNC_symbol'], include='all')
        line['raw'] = entry

        yield line

def main(argv):

    header = [ 'Disease_trivial_name', 'HGNC_symbol', 'OMIM_morbid', 'Gene_locus', 'Chromosome', 'Clinical_db_gene_annotation' ]

    infile = open(argv[0])

    raw_data = ( line.strip() for line in infile) # sluuuurp file
    next(raw_data) # remove the header
    parsable_data = ( line.split("\t") for line in raw_data )


    data = list2dict(header, parsable_data)
    json_data = get_raw_omim(data)

    all_json = [ line['raw'] for line in json_data ]

    print('[ {} ]'.format(','.join(all_json)))

if __name__ == '__main__':
    main(sys.argv[1:])
