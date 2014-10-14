#!/usr/bin/env python
# encoding: utf-8

# Usage:
#   create_list.py <csv-infile> [--hgnc] > gene_list.csv
#
# Arguments:
#   csv-infile a one column csv file with either HGNC or EnsEMBL IDs
#
# Options:
#   --hgnc: csv-infile is filled with HGNC identifiers
#
# Will output a genelist with the right right amount of column filled in with #NA values.
# The right column with identifiers from the csv-infile will be filled in.
from __future__ import print_function
import sys
import argparse


def main(argv):
    """
    :argv: input arguments of the script
    :returns: pass

    """

    # set up the argparser
    parser = argparse.ArgumentParser(description='Will output a genelist with \
            the right right amount of column filled in with #NA values. The \
            right column with identifiers from the csv-infile will be filled in.')
    parser.add_argument('--hgnc', dest='hgnc', default=False, action='store_true', help='Indicate that the csv-infile is filled with HGNC identifiers')
    parser.add_argument('infile', type=argparse.FileType('r'), help='a one column csv file with either HGNC or EnsEMBL IDs')
    args = parser.parse_args(argv)

    infile = args.infile
    identifiers = [ line.strip() for line in infile.readlines() ] # trim line
    identifiers = [ line for line in identifiers if line ]        # skip blank lines

    # print header
    print("#Chromosome\tGene_start\tGene_stop\tHGNC_ID\tDisease_group_pathway\tProtein_name\tSymptoms\tBiochemistry\tImaging\tDisease_trivial_name\tTrivial_name_short\tGenetic_model\tOMIM_gene\tOMIM_morbid\tGene_locus\tGenome_build\tUniPort_ID\tEnsembl_gene_id\tEnsemble_transcript_ID\tRed_pen\tDatabase")
    if args.hgnc:
        for identifier in identifiers:
            print("#NA\t" * 3, identifier, "\t#NA" * 17)
    else:
        for identifier in identifiers:
            print("#NA\t" * 17, identifier, "\t#NA" * 3)

if __name__ == '__main__':
    main(sys.argv[1:])
