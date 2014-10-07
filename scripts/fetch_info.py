#!/usr/bin/env python
# encoding: utf-8

import sys
import pymysql
import argparse
import re

gl_header=['#Chromosome', 'Gene_start', 'Gene_stop', 'HGNC_ID', 'Disease_group_pathway', 'Protein_name', 'Symptoms', 'Biochemistry', 'Imaging', 'Disease_trivial_name', 'Trivial_name_short', 'Genetic_model', 'OMIM_gene', 'OMIM_morbid', 'Gene_locus', 'Genome_build', 'UniPort_ID', 'Ensembl_gene_id', 'Ensemble_transcript_ID', 'Red_pen', 'Database']

HGNC_FIELD = 3
ENSEMBL_FIELD = 17

def print_header(header=gl_header):
    print_line(header)

def print_line(line):
    print("\t".join(line))

def query_hgnc(data):
    conn = pymysql.connect(host='ensembldb.ensembl.org', port=5306, user='anonymous', db='homo_sapiens_core_75_37')
    cur = conn.cursor()
    for idx, line in enumerate(data):
        print("Getting '%s' ... " % line[HGNC_FIELD])
        cur.execute("select seq_region.name AS Chromosome, g.seq_region_start, g.seq_region_end, x.display_label, '#NA', '#NA',  '#NA',  '#NA',  '#NA',  '#NA',  '#NA', '#NA', '#NA', '#NA', '#NA', 'GRCh38', '#NA', g.stable_id, '#NA','#NA',  'Database' from gene g join xref x on x.xref_id = g.display_xref_id join seq_region using (seq_region_id) where x.display_label = %s", line[HGNC_FIELD])
    
        for r in cur.fetchall():
            print(r)
    return data

def query_ens(data):
    pass



def main(argv):
    # set up the argparser
    parser = argparse.ArgumentParser(description='Queries EnsEMBL and fills in the blanks of a gene list. Only columns headers found in gl_headers will be used')
    parser.add_argument('--hgnc', dest='hgnc', default=False, action='store_true', help='Indicate that the csv-infile''s main identifiers are HGNC identifiers')
    parser.add_argument('infile', type=argparse.FileType('r'), help='the csv file')
    args = parser.parse_args(argv)

    # read in the CSV file
    csvfile = args.infile
    data = [ line.strip() for line in csvfile.readlines() ] # sluuuurp
    data = [ line.split("\t") for line in data ]
    
    # clean up the input
    for idx, line in enumerate(data):
        # replace white space seperated comma's with just a comma
        # replace ; with comma
        # remove leading and trailing white space
        line = [ re.sub(r'\s*[,;]\s*', ',', column.strip()) for column in line ] 
        data[idx] = line
    header = data.pop(0)

    if args.hgnc:
        data = query_hgnc(data)
    else:
        data = query_ens(data)

    print_header()
    for line in data:
        print_line(line)

if __name__ == '__main__':
    main(sys.argv[1:])
