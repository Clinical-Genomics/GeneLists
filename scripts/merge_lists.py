#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function
import sys
import argparse
import subprocess
import os

def getgittag(filename):
    """Gets the current version of a gene list

    Args:
        filename (str): the name of the gene list

    Returns (str): a version (tag) of the gene list

    """
    cwd = os.getcwd()
    os.chdir(os.path.dirname(filename))
    tag = subprocess.check_output(['git', 'describe']).decode('utf-8').strip()
    os.chdir(cwd)

    return tag


def main(argv):
    parser = argparse.ArgumentParser(description='Merge gene lists. Will only output HGNC_ID, EnsEMBL_gene_id and Database columns.')
    parser.add_argument('infiles', nargs="+", type=argparse.FileType('r'), help='')
    parser.add_argument('-d', '--database', nargs='*', dest='database', action='append', help='only take HGNC_IDs from this database.')
    args = parser.parse_args(argv)

    databases = []
    if len(args.database):
        databases = [ db[0] for db in args.database ]

    versions = {} # Filename => { Database => Version }
    data = {} # HGNC_ID => {'HGNC_ID' => '', 'EnsEMBLid' => [], 'Databases' => () }
    for infile in args.infiles:
        versions[infile.name] = {}
        for line in infile:
            if line.startswith('#'): continue
            line = line.split("\t")
            hgnc_id = line[3].strip()
            ensEMBLid = line[17].strip()
            red_pen = line[19].strip()
            line_databases = line[20].strip().split(',')
            dis_ass_trans = ''
            if len(line) == 22:
                dis_ass_trans = line[21].strip()

            # skip if we are whitelisting dbs
            if len(databases):
                set_databases = set(line_databases).intersection(databases)
                if len(set_databases):
                    line_databases = list(set_databases)
                else:
                    continue

            # init
            if hgnc_id not in data:
                data[hgnc_id] = {'HGNC_ID': '', 'EnsEMBLid': [], 'Databases': [], 'red_pen': '', 'dis_ass_trans': '' }

            # fill
            data[hgnc_id]['HGNC_ID'] = hgnc_id
            data[hgnc_id]['EnsEMBLid'].append(ensEMBLid)
            data[hgnc_id]['Databases'] += line_databases
            data[hgnc_id]['red_pen'] = red_pen
            data[hgnc_id]['dis_ass_trans'] = dis_ass_trans

            # fill versions dict
            for database in line_databases:
                if database not in versions[infile.name]:
                    version = getgittag(infile.name)
                    versions[infile.name][database] = getgittag(infile.name)

    for filename, database_version in versions.items():
        for database, version in database_version.items():
            print('##Database=<ID=%s,Version=%s,Acronym=%s,Clinical_db_genome_build=GRCh37.p13' % (os.path.basename(filename), version, database))

    print('HGNC_ID	Ensembl_gene_id	Clinical_db_gene_annotation	Reduced_penetrance	Disease_associated_transcript')
    for line in data.values():
        line['Databases'].append('FullList')
        line_dbs = ','.join(list(set(line['Databases'])))
        if line['dis_ass_trans'] and line['dis_ass_trans'] != '#NA':
            dis_ass_trans = '%s:%s|' % (line['HGNC_ID'], line['dis_ass_trans'])
        else:
            dis_ass_trans = line['dis_ass_trans']
        for ensEMBLid in set(line['EnsEMBLid']):
            print('%s\t%s\t%s\t%s\t%s' % (line['HGNC_ID'], ensEMBLid, line_dbs, line['red_pen'], dis_ass_trans))

if __name__ == '__main__':
    main(sys.argv[1:])
