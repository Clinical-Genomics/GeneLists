""" Merge the panels of gene lists into one list """
# encoding: utf-8

from __future__ import print_function
import os
import datetime

from ..utils.git import getgitlastmoddate, getgittag
from ..utils.acronyms import Acronyms

def merge_panels(infiles, databases):
    """ Merge the panels of gene lists into one list. """

    acronyms = Acronyms(os.path.dirname(os.path.dirname(os.path.realpath(infiles[0].name))))

    versions = {} # Filename => { Database => { 'Version': Version, 'Date': Date } }
    data = {} # HGNC_symbol => {'HGNC_symbol' => '', 'EnsEMBLid' => [], 'Databases' => () }

    columns = ['Chromosome', 'HGNC_symbol', 'Ensembl_gene_id', 'Clinical_db_gene_annotation',
               'Reduced_penetrance', 'Disease_associated_transcript', 'Phenotypic_disease_model',
               'Genetic_disease_model', 'OMIM_morbid', 'Database_entry_version', 'Curator',
               'Alias', 'Group_or_Pathway', 'Mosaicism', 'Comments']

    list_columns = ['Genetic_disease_model', 'Ensembl_gene_id', 'Clinical_db_gene_annotation']

    for infile in infiles:
        versions[infile.name] = {}
        for line in infile:
            # deal with headers
            if line.startswith('##'):
                continue # skip comments
            if line.startswith('#Chromosome'): # get the header
                line = line.lstrip('#')
                header = line.split("\t")
                continue

            line = line.split("\t")
            line = dict(zip(header, line)) # get me a nice mapping

            # sanitize all columns
            for column in columns:
                if column not in line:
                    line[column] = ''
                line[column] = line[column].strip()
                if ':' in line[column]:
                    line[column] = line[column].split(':')[1]

            # the models can be multiple, so make it into a list
            for column in list_columns:
                line[column] = line[column].split(',')

            # skip if we are whitelisting dbs
            if len(databases):
                set_databases = set(line['Clinical_db_gene_annotation']).intersection(databases)
                if len(set_databases):
                    line['Clinical_db_gene_annotation'] = list(set_databases)
                else:
                    continue

            hgnc_id = line['HGNC_symbol']

            # init
            if hgnc_id not in data:
                data[hgnc_id] = dict(zip(columns, ['' for i in range(len(columns))]))
                for column in list_columns:
                    data[hgnc_id][column] = []

            # fill
            for column in columns:
                if column in list_columns:
                    data[hgnc_id][column].append(line[column])
                data[hgnc_id][column] = line[column]

            # fill versions dict
            for database in line['Clinical_db_gene_annotation']:
                if database == 'OMIM':
                    if not infile.name.endswith('OMIM.txt'):
                        continue
                if database not in versions[infile.name]:
                    full_mod_date = getgitlastmoddate(infile.name, '%c')
                    if not full_mod_date: # ok, we haven't saved this list yet
                        mod_date = datetime.datetime.now().strftime('%Y%m%d')
                    else:
                        mod_date = datetime.datetime.\
                                   strptime(full_mod_date, '%c').\
                                   strftime('%Y%m%d')
                    version = getgittag(infile.name, full_mod_date)
                    full_name = acronyms[database]
                    versions[infile.name][database] = {'Version': version, 'Date': mod_date,
                                                       'Fullname': full_name}

    for filename, database_version in versions.items():
        for database, version_date in database_version.items():
            print('##Database=<ID=%s,Version=%s,Date=%s,Acronym=%s,Complete_name=%s,Clinical_db_genome_build=GRCh37.p13' % (os.path.basename(filename), version_date['Version'], version_date['Date'], database, version_date['Fullname']))

    print('\t'.join(columns))
    for line in data.values():
        if len(databases) > 2:
            line['Clinical_db_gene_annotation'].append('FullList')

        for column in list_columns:
            line[column] = ','.join(list(set(line[column])))

        if line['Disease_associated_transcript'] == 'unknown':
            line['Disease_associated_transcript'] = ''
        for ensEMBLid in line['Ensembl_gene_id'].split(','):
            print('\t'.join([ line[column] for column in columns ]))
