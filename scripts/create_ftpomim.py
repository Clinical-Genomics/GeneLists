#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function
import sys
import argparse
import re

omim_header=['Disease_trivial_name', 'HGNC_symbol', 'OMIM_morbid', 'Gene_locus', 'Chromosome', 'Clinical_db_gene_annotation']

def print_header(header=omim_header):
    """prints column headers
    Args:
        header (list, optional): a list of strings
    Returns:
        pass
    Note:
        Prints to STDOUT
    """
    print("\t".join(header))

def print_line(line):
    """Prints a line based on the order of the gl_headers

    Args:
        line (list): list with values for one line of a gene list.

    Returns:
        pass

    Note:
        Will print the STDOUT
    """
    print("\t".join( ( str(part) for part in line) ))

def pick_hgnc_symbol(data):
    """Removes all but one HGNC symbol. It will keep the first one.

    Args:
        data (list): list with values for one line of a gene list.

    Yields:
        list: with all but one hgnc symbol removed.
    """
    for line in data:
        hgnc_symbol = line[1].split(',')
        hgnc_symbol = hgnc_symbol[0].strip() if hgnc_symbol else None
        line[1] = hgnc_symbol
        yield line

def remove_duplicates(data):
    """Removes lines with the same OMIM ID. First line encountered is kept.

    Args:
        data (list): list with values for one line of a gene list.

    Yields:
        list: only first line with the OMIM ID
    """
    omim_ids = {} # omim_id => 1
    for line in data:
        omim_id = line[2]
        if omim_id not in omim_ids:
            yield line
        omim_ids[ omim_id ] = 1

def add_chromosome_number(data):
    """TODO: Docstring for add_chromosome_number.

    Args:
        data (list): list with values for one line of a gene list.

    Yields:
        list: with one extra column for the chromosome number
    """
    for line in data:
        locus=line[3]
        chromosome=''
        if any([x for x in locus if x in ('q', 'p')]): # check if the locus has an q or p
            chromosome = re.compile('p|q').split(locus)[0]
        elif locus.startswith('Chr.'):
            chromosome = locus.replace('Chr.', '')

        line.append(chromosome)
        yield line

def add_database(data, database='OMIM'):
    """TODO: Docstring for add_database.

    Args:
        data (list): list with values for one line of a gene list.

    Yields:
        list: with one extra column for the database
    """

    for line in data:
        line.append(database)
        yield line

def remove_unconfirmed(data):
    """TODO: Docstring for remove_.

    Args:
        arg1 (TODO): TODO

    Returns: TODO

    """
    for line in data:
        if line[0].startswith('?'): continue
        if line[0].startswith('{'): continue
        yield line

def main(argv):
    # set up the argparser
    parser = argparse.ArgumentParser(description='Creates a gene list from an OMIM morbid file')
    parser.add_argument('infile', type=argparse.FileType('r'), help='the OMIM morbid file with columns: phenotype, HGNC symbols, OMIM ID, gene locus')
    parser.add_argument('-d', default='OMIM', dest='database', help='The value for the database field. Defaults to "OMIM"')
    args = parser.parse_args(argv)

    omimfile = args.infile
    database = args.database
    raw_data = ( line.strip() for line in omimfile) # sluuuurp
    parsable_data = ( line.split("|") for line in raw_data )

    print_header()
    [ print_line(line) for line in \
      add_database(
      add_chromosome_number(
      remove_unconfirmed(
      remove_duplicates(
        parsable_data
      ))), database)
    ]
    #hgnc_data = pick_hgnc_symbol(uniq_data)

if __name__ == '__main__':
    main(sys.argv[1:])
