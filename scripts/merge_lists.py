#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function
import sys
import argparse
import subprocess
import os

from .git import getgitlastmoddate, getgittag
from .acronyms import Acronyms

def main(argv):
    parser = argparse.ArgumentParser(description='Merge gene lists. Will only output HGNC_symbol, EnsEMBL_gene_id and Database columns.')
    parser.add_argument('infiles', nargs="+", type=argparse.FileType('r'), help='')
    parser.add_argument('-d', '--database', nargs='*', dest='database', action='append', help='only take HGNC_symbols from this database.')
    args = parser.parse_args(argv)

    databases = []
    if len(args.database):
        databases = [ db[0] for db in args.database ]

    acronyms = Acronyms(os.path.dirname(os.path.dirname(os.path.realpath(args.infiles[0].name))))

    versions = {} # Filename => { Database => { 'Version': Version, 'Date': Date } }
    data = {} # HGNC_symbol => {'HGNC_symbol' => '', 'EnsEMBLid' => [], 'Databases' => () }
    for infile in args.infiles:
        versions[infile.name] = {}
        for line in infile:
            if line.startswith('#'): continue # skip comments
            line = line.split("\t")
            if len(line) < 21: continue # crude way of omitting non-gene lists
            hgnc_id = line[3].strip()
            phenotypic_disease_model = line[10].strip() # keep manual annotations
            ensEMBLid = line[14].strip()
            red_pen = line[16].strip()
            line_databases = line[17].strip().split(',')
            dis_ass_trans = line[18].strip() 

            # skip if we are whitelisting dbs
            if len(databases):
                set_databases = set(line_databases).intersection(databases)
                if len(set_databases):
                    line_databases = list(set_databases)
                else:
                    continue

            # init
            if hgnc_id not in data:
                data[hgnc_id] = {'HGNC_symbol': '', 'EnsEMBLid': [], 'Databases': [], 'red_pen': '', 'dis_ass_trans': '' }

            # fill
            data[hgnc_id]['HGNC_symbol'] = hgnc_id
            data[hgnc_id]['EnsEMBLid'].append(ensEMBLid)
            data[hgnc_id]['Databases'] += line_databases
            data[hgnc_id]['red_pen'] = red_pen
            data[hgnc_id]['dis_ass_trans'] = dis_ass_trans
            data[hgnc_id]['Phenotypic_disease_model'] = phenotypic_disease_model

            # fill versions dict
            for database in line_databases:
                if database not in versions[infile.name]:
                    version = getgittag(infile.name)
                    mod_date = getgitlastmoddate(infile.name)
                    full_name = acronyms[database]
                    versions[infile.name][database] = { 'Version': version, 'Date': mod_date, 'Fullname': full_name }

    for filename, database_version in versions.items():
        for database, version_date in database_version.items():
            print('##Database=<ID=%s,Version=%s,Date=%s,Acronym=%s,Complete_name=%s,Clinical_db_genome_build=GRCh37.p13' % (os.path.basename(filename), version_date['Version'], version_date['Date'], database, version_date['Fullname']))

    print('HGNC_symbol	Ensembl_gene_id	Clinical_db_gene_annotation	Reduced_penetrance	Disease_associated_transcript	Phenotypic_disease_model')
    for line in data.values():
        line['Databases'].append('FullList')
        line_dbs = ','.join(list(set(line['Databases'])))
        if line['dis_ass_trans'] == 'unknown': line['dis_ass_trans'] = '#NA'
        if line['dis_ass_trans'] and line['dis_ass_trans'] != '#NA':
            dis_ass_trans = '%s:%s' % (line['HGNC_symbol'], line['dis_ass_trans'])
        else:
            dis_ass_trans = line['dis_ass_trans']
        for ensEMBLid in set(line['EnsEMBLid']):
            print('%s\t%s\t%s\t%s\t%s\t%s' % (line['HGNC_symbol'], ensEMBLid, line_dbs, line['red_pen'], dis_ass_trans, line['Phenotypic_disease_model']))

if __name__ == '__main__':
    main(sys.argv[1:])
