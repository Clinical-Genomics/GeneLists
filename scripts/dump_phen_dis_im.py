#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function
import sys

from .omim import OMIM

def p(line):
    """print only if the verbose switch has been set

    Args:
            line (str): line to print to STDOUT

    Returns:
            pass
    """
    print('\033[93m', '>>> ', line, '\033[0m')


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

def yield_genes(data):
    """
    Only yield gene and gene/phenotypes

    Args:
        data (list of dicts): Inner dict represents a row in a gene list
    Yields:
        dict: with 
    """
    for line in data:
        if line['Type'] in ('gene', 'gene/phenotype'):
            yield line

def query_omim(data):
    """Queries OMIM to fill in the inheritance models

    Args:
            data (list of dicts): Inner dict represents a row in a gene list

    Yields:
            dict: with Inheritance_models and 'Phen MIM number'
    """

    omim = OMIM(api_key='3BA13F7F94F69F7F1B1AC3FCB896ECB6FA8B8D42')
    for line in data:
        if 'Approved Gene Symbol' in line:
            entry = omim.gene(line['Approved Gene Symbol'])

            phenotypic_disease_models = omim.parse_phenotypic_disease_models(entry['phenotypes'])

            # extract the inheritance model
            line_phenotypic_disease_models = []

            # if any inheritance models and omim numbers are present, use them!
            for phen_mim_number, inheritance_models_descr in phenotypic_disease_models.items():
                if phen_mim_number is not None:
                    inheritance_models_str = ''

                    if inheritance_models_descr['models'] is not None:
                        inheritance_models_str = ','.join(inheritance_models_descr['models'])
                        line['Inheritance models'] = inheritance_models_str
                        line['Phen MIM number'] = phen_mim_number
                        line['Description'] = inheritance_models_descr['descr']
                        if (str(line['MIM Number']) != str(entry['mim_number'])):
                            p('%s > %s client OMIM number differs from OMIM query' % (entry['mim_number'], line['MIM Number']))
                        yield line

def main(argv):

    header = ['MIM Number', 'Type', 'Entrez Gene ID', 'Approved Gene Symbol', 'Ensembl Gene ID']

    infile = open(argv[0]) # needs to mim2gene file

    raw_data = (line.strip() for line in infile) # sluuuurp file
    next(raw_data) # remove the header
    parsable_data = (line.split("\t") for line in raw_data)

    data = list2dict(header, parsable_data)
    lines = query_omim(yield_genes(data))

    for line in lines:
        print('{}\t{}\t{}\t{}'.format(line['MIM Number'], line['Phen MIM number'], line['Inheritance models'], line['Description']))

if __name__ == '__main__':
    main(sys.argv[1:])
