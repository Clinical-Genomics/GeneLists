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
            if len(line) < 21:
                continue # crude way of omitting non-gene lists
            line = dict(zip(header, line)) # get me a nice mapping

            hgnc_id = line['HGNC_symbol'].strip()
            chromosome = line['Chromosome'].strip()
            # keep manual annotations
            phenotypic_disease_model = line['Phenotypic_disease_model'].strip()
            ensEMBLid = line['Ensembl_gene_id'].strip()
            red_pen = line['Reduced_penetrance'].strip()
            line_databases = line['Clinical_db_gene_annotation'].strip().split(',')
            dis_ass_trans = line['Disease_associated_transcript'].strip()
            genetic_disease_model = line['Genetic_disease_model'].strip() \
                                    if 'Genetic_disease_model' in line else ''
            omim_morbid = line['OMIM_morbid'].strip()

            # sanitize omim morbid
            if ':' in omim_morbid:
                omim_morbid = omim_morbid.split(':')[1]

            # skip if we are whitelisting dbs
            if len(databases):
                set_databases = set(line_databases).intersection(databases)
                if len(set_databases):
                    line_databases = list(set_databases)
                else:
                    continue

            # init
            if hgnc_id not in data:
                data[hgnc_id] = {'HGNC_symbol': '', 'EnsEMBLid': [], 'Databases': [],
                                 'red_pen': '', 'dis_ass_trans': '', 'genetic_disease_model': '' }

            # fill
            data[hgnc_id]['HGNC_symbol'] = hgnc_id
            data[hgnc_id]['Chromosome'] = chromosome
            data[hgnc_id]['EnsEMBLid'].append(ensEMBLid)
            data[hgnc_id]['Databases'] += line_databases
            data[hgnc_id]['red_pen'] = red_pen
            data[hgnc_id]['dis_ass_trans'] = dis_ass_trans
            data[hgnc_id]['Phenotypic_disease_model'] = phenotypic_disease_model
            data[hgnc_id]['genetic_disease_model'] = genetic_disease_model
            data[hgnc_id]['OMIM_morbid'] = omim_morbid

            # fill versions dict
            for database in line_databases:
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

    print('Chromosome	HGNC_symbol	Ensembl_gene_id	Clinical_db_gene_annotation	Reduced_penetrance	Disease_associated_transcript	Phenotypic_disease_model	Genetic_disease_model	OMIM_morbid')
    for line in data.values():
        if len(databases) > 2:
            line['Databases'].append('FullList')
        line_dbs = ','.join(list(set(line['Databases'])))
        if line['dis_ass_trans'] == 'unknown':
            line['dis_ass_trans'] = ''
        if line['dis_ass_trans'] and line['dis_ass_trans'] != '#NA':
            dis_ass_trans = '%s:%s' % (line['HGNC_symbol'], line['dis_ass_trans'])
        else:
            dis_ass_trans = line['dis_ass_trans']
        for ensEMBLid in set(line['EnsEMBLid']):
            print('\t'.join((line['Chromosome'], line['HGNC_symbol'], ensEMBLid, line_dbs,
                             line['red_pen'], dis_ass_trans, line['Phenotypic_disease_model'],
                             line['genetic_disease_model'], line['OMIM_morbid'])))
