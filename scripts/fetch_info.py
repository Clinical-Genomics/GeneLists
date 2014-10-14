#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function
import sys
import pymysql
import argparse
import re

gl_header=['Chromosome', 'Gene_start', 'Gene_stop', 'HGNC_ID', 'Disease_group_pathway', 'Protein_name', 'Symptoms', 'Biochemistry', 'Imaging', 'Disease_trivial_name', 'Trivial_name_short', 'Genetic_model', 'OMIM_gene', 'OMIM_morbid', 'Gene_locus', 'Genome_build', 'UniPort_ID', 'Ensembl_gene_id', 'Ensemble_transcript_ID', 'Red_pen', 'Database']

def print_header(header=gl_header):
    """Prints the gl_header

    Args:
        header (list, optional): a list of strings
    Returns:
        pass
    Note:
        Prints to STDOUT
    """
    print('#' + "\t".join(gl_header))

def print_line(line):
    """Prints a line based on the order of the gl_headers

    Args:
        line (dict): dict with values for one line of a gene list. All values of gl_headers should be present as keys

    Returns:
        pass

    Note:
        Will print the STDOUT
    """
    ordered_line = list()
    for column_name in gl_header:
        ordered_line.append(str(line[column_name]))
    print("\t".join(ordered_line))

def query(data, keys):
    """Queries EnsEMBL. Parameters are HGNC_ID and/or Ensembl_gene_id, whatever is available. Data from EnsEMBLdb will overwrite the client data.
    A(n) identifier(s) should yield one result from EnsEMBLdb. It will be reported if a(n) identifier(s) don't yield any or multiple results.
    
    Args:
        data (list of dicts): representing the lines and columns in a gene list. The keys of the dicts must match the column names of the EnsEMBLdb query.
        keys (list): A list of available parameters: HGNC_ID and/or Ensembl_gene_id.
            TODO: could this now be autodetected?
    Yields:
        dict: a row with data from ensEMBLdb filled in.
    
    """
    conn = pymysql.connect(host='ensembldb.ensembl.org', port=5306, user='anonymous', db='homo_sapiens_core_75_37')
    cur = conn.cursor(pymysql.cursors.DictCursor)

    base_query = "select g.seq_region_start AS Gene_start, g.seq_region_end AS Gene_stop, x.display_label AS HGNC_ID, g.stable_id AS Ensembl_gene_id, seq_region.name AS Chromosome from gene g join xref x on x.xref_id = g.display_xref_id join seq_region using (seq_region_id)"
    keys_conds = { 'HGNC_ID': 'x.display_label', 'Ensembl_gene_id': 'g.stable_id' }

    for line in data:
        conds = [ "%s = %%s" % keys_conds[ key ] for key in keys if key in line ]
        cond_values = [ line[ key ] for key in keys if key in line ]

        query = "%s where %s" % ( base_query, " and ".join(conds) )
        cur.execute(query, cond_values)

        rs = cur.fetchall() # result set

        # O-oh .. for now this still means manual intervention!
        if len(rs) == 0:
            print("Getting '%s' ... " % cond_values)
            print('Not found!')
        elif len(rs) > 1:
            print("Getting '%s' ... " % cond_values)
            print('Multiple entries found!')
        else:
            for entry in rs:
                yield merge_line(entry, line)

def merge_line(ens, client):
    """Will merge dict ens (EnsEMBL data) with client (data). ens will take precedence over client. Changes will be reported.

    Args:
        ens (dict): dict with values for one line of a gene list
        client (dict): with values for one line of a gene list
    Yields:
        dict: merged ens and client dict

    """
    for key, value in client.items():
        if key in ens:
            if str(ens[key]) != str(value):
                print("%s > %s: ens '%s' diff from client '%s'" % (ens['Ensembl_gene_id'], key, ens[key], value))
        #    else:
        #        print("%s: ens '%s' eq to client '%s'" % (key, ens[key], value))
        #else:
        #    print("%s not in ens!" % key)

    merged = client.copy()
    merged.update(ens)
    return merged

def fill(data):
    """Fills in the blanks with '#NA'

    Args:
        data (list of dicts): representing the lines and columns in a gene list
    Yields:
        dict: with all missing columns filled in with #NA

    """
    defaults=dict((column_name, '#NA') for column_name in gl_header)
    for line in data:
        d = defaults.copy()
        d.update(line)
        yield d

def munge(data):
    """Make sure the data we got from EnsEMBL is good enough for the gene lists
        # swap start and stop gene coordinates if start if bigger than stop

    Args:
        data (list of dicts): representing the lines and columns in a gene list
    Yields:
        dict: a row with some munged data

    """
    # swap coordinates if start > stop
    for line in data:
        if line['Gene_start'] > line['Gene_stop']:
            line['Gene_stop'], line['Gene_start'] = line['Gene_start'], line['Gene_stop']
        yield line

def cleanup(data):
    """Will clean the data according to rules set by MIP
        # replace white space seperated comma's with just a comma
        # replace ; with comma
        # remove leading and trailing white space

    Args:
        data (list of lists): Inner list represents a row in a gene list
    Yield:
        list: a cleaned up row

    """
    for line in data:
        yield [ re.sub(r'\s*[,;]\s*', ',', column.strip()) for column in line ]

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

def main(argv):
    # set up the argparser
    parser = argparse.ArgumentParser(description='Queries EnsEMBL and fills in the blanks of a gene list. Only columns headers found in gl_headers will be used')
    parser.add_argument('infile', type=argparse.FileType('r'), help='the tsv file with correct headers')
    args = parser.parse_args(argv)

    # read in the TSV file
    tsvfile = args.infile
    raw_data = ( line.strip() for line in tsvfile ) # sluuuurp
    parsable_data = ( line.split("\t") for line in raw_data )

    # clean up the input
    clean_data = cleanup(parsable_data)

    # list to dict
    header = next(clean_data) # get the header
    dict_data = list2dict(header, clean_data)

    # fill in missing blanks
    ensembld_data = query(dict_data, ['HGNC_ID', 'Ensembl_gene_id'])

    # clean up the data from EnsEMBL a bit
    munged_data = munge(ensembld_data)

    # fill in missing values with #NA
    completed_data = fill(munged_data)

    # print the gene list
    print_header()
    for line in completed_data:
        print_line(line)

if __name__ == '__main__':
    main(sys.argv[1:])
