#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function
import sys
import argparse
import re

gl_header=['Chromosome', 'Gene_start', 'Gene_stop', 'HGNC_symbol', 'Protein_name', 'Symptoms', 'Biochemistry', 'Imaging', 'Disease_trivial_name', 'Trivial_name_short', 'Phenotypic_disease_model', 'OMIM_morbid', 'Gene_locus', 'UniProt_id', 'Ensembl_gene_id', 'Ensemble_transcript_ID', 'Reduced_penetrance', 'Clinical_db_gene_annotation', 'Disease_associated_transcript', 'Ensembl_transcript_to_refseq_transcript', 'Gene_description']
mandatory_fields = {
    'Clinical_db_gene_annotation': re.compile(r'.+'),
    'Chromosome': re.compile(r'([\dXY]|MT)+'),
    'Gene_start': re.compile(r'\d+'),
    'Gene_stop': re.compile(r'\d+'),
    'HGNC_symbol': re.compile(r'.+'),
#    'Genetic_disease_model': re.compile(r'.+'),
#    'Gene_locus': re.compile(r'.+'),
    'Ensembl_gene_id': re.compile(r'ENSG\d{11}'),
}

# holds regex's with forbidden chars per column
forbidden_chars = {
    # Can't have empty inheritance models:
    # Fobidden: SYMBOL:>, SYMBOL:OMIM||, SYMBOL:OMIM>, SYMBOL:OMIM>AR|>
    'Phenotypic_disease_model': re.compile(r'(:>|>$|:\d+>\||\|>)')
}

line_nr = 0 # hold on to the line nr
warned = 0 # exit code

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

def warn(msg):
    """Pretty print a warning message

    Args:
        msg (str): The message to print

    """
    global line_nr
    global warned
    print("#{}: {}".format(line_nr, msg))
    warned = 1

def inc_line_nr(lines):
    """Increments the global line_nr for each passing line

    Args:
        lines (list of dicts): each dict contains a dict with gl_header as keys
    Returns: pass

    """
    global line_nr
    for line in lines:
        line_nr += 1
        yield line

def check_forbidden_chars(lines):
    """Checks the precense of forbidden chars

    Args:
        lines (list of dicts): each dict contains a dict with gl_header as keys

    yields: a line of the lines

    """
    for line in lines:
        for field, forbidden_re in forbidden_chars.items():
            if re.search(forbidden_re, line[field]):
                warn("'{}' ('{}') has a forbidden char combination '{}'".format(field, line[field], forbidden_re))
        yield line

def check_delimiter(lines, delimiters=re.compile(r', '), whitelist_keys=[]):
    """Check for the presence of a delimiter. Will use the 'warn' function to print presence of a delimiter

    Args:
        lines (list of dicts): each dict contains a dict with gl_header as keys
        delimiters (re.compile): a compiled regex with what to match

    yields: a line of the lines

    """
    for line in lines:
        for key, value in line.items():
            if key in whitelist_keys: continue
            if re.search(delimiters, value):
                warn("{} ('{}') holds a forbidden delimiter".format(key, value))
        yield line

def check_trimming(lines):
    """Check if fields are trimmed

    Args:
        lines (list of dicts): each dict contains a dict with gl_header as keys

    yields: a line of the lines
    """
    for line in lines:
        for key, value in line.items():
            if value != value.strip():
                warn("{} ('{}') is not trimmed!".format(key, value))
        yield line

def check_nr_fields(lines):
    """Checks if the nr of fields == nr of headers

    Args:
        lines (list of dicts): each dict contains a dict with gl_header as keys

    yields: a line of the lines
    """
    global gl_header
    len_header = len(gl_header)
    for line in lines:
        len_line = len(line)
        if len_line != len_header:
            warn("Len fields ({}) != len headers ({})".format(len_line, len_header))
        yield line

def check_mandatory_fields(lines):
    """Check if certain fields are filled in

    Args:
        lines (list of dicts): each dict contains a dict with gl_header as keys

    yields: a line of the lines
    """
    for line in lines:
        for mandatory_field, mandatory_re in mandatory_fields.items():
            if not re.search(mandatory_re, line[mandatory_field]):
                warn("'{}' ('{}') fails mandatory requirement '{}'".format(mandatory_field, line[mandatory_field], mandatory_re))
        yield line

def check_coordinates(lines):
    """Checks if gene length is above zero

    Args:
        lines (list of dicts): each dict contains a dict with gl_header as keys

    yields: a line of the lines
    """
    for line in lines:
        if int(line['Gene_stop']) - int(line['Gene_start']) <= 0:
            warn('Gene coordinates are not above zero.')
        yield line

def check_duplicates(lines):
    """Check for duplicates based on HGNC_symbol, EnsEMBL_gene_id and coordinates

    Args:
        lines (list of dicts): each dict contains a dict with gl_header as keys

    yields: a line of the lines
    """
    fields = {
        'HGNC_symbol': {}, # HGNC_symbol: line nr
        'Ensembl_gene_id': {},
    }
    Gene_start = {}
    Gene_stop = {}
    for line in lines:
        for field_key in fields.keys():
            if line[field_key] in fields[field_key]:
                warn("'{}' already listed at #{}".format(line[field_key], fields[field_key][ line[field_key] ]))
            fields[field_key][ line[field_key] ] = line_nr

        if line['Gene_start'] in Gene_start and line['Gene_stop'] in Gene_stop:
            warn("'{}-{}' already listed at #{}".format(line['Gene_start'], line['Gene_stop'], Gene_start[ line['Gene_start'] ]))
        Gene_start[ line['Gene_start'] ] = line_nr
        Gene_stop[  line['Gene_stop']  ] = line_nr

        yield line

def main(argv):
    """Main program

    Args:
        argv (list): input params

    Returns: 0 on success, otherwise error code

    """
    parser = argparse.ArgumentParser(description='Sanity checks a gene list')
    parser.add_argument('infile', type=argparse.FileType('r'), help='the tsv file with correct headers')
    args = parser.parse_args(argv)

    infile = args.infile

    lines = (line.strip('\r\n') for line in infile) # sluuuurp

    # get the line in parts
    parsable_data = ( line.split('\t') for line in lines )

    # skip parsing of leading comments
    comments = []
    line = next(parsable_data)
    global line_nr
    line_nr = 1
    while line[0].startswith('##'):
        comments.append(line)
        line = next(parsable_data)
        line_nr += 1

    # OK - start the sanity checks

    # header starts with #
    header = line # get the header
    if not header[0].startswith('#'):
        warn('Header does not start with #')
    else:
        header[0] = header[0][1:] # remove the leading #

    dict_data = list2dict(header, parsable_data)

    # are all mandatory headers there?
    global gl_header
    global mandatory_fields
    header_diff = set(mandatory_fields.keys()).difference(header)
    if len(header_diff) > 0: warn("Missing columns {}".format(header_diff))

    # header should not contain white space
    for head in header:
        if re.search(' ', head) != None:
            warn("Header '{}' contains white space".format(head))

    # file contains ';' delimiter
    # elements in a field should be seperated by ',', not ', '
    # fields should be trimmed
    # #fields = #fields in header
    # following fields should be filled in:
    # - clinical_db_gene_annotation
    # - clinical_db_genome_build
    # - chromosome
    # - gene_start
    # - gene_stop
    # - HGNC_symbol
    # - Genetic_model
    # - Gene_locus
    # - Ensembl_gene_id
    # length of coordinates should be above zero
    # no duplicates based on coordinates tuple, EnsEMBL gene ID or HGNC symbol
    [ line for line in \
        check_duplicates(
        check_coordinates(
        check_mandatory_fields(
        check_forbidden_chars(
        check_nr_fields(
        check_trimming(
        check_delimiter(
        inc_line_nr(
            dict_data
        ))))))))
    ]

    return warned

if __name__ == '__main__':
    main(sys.argv[1:])
