#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function
import os
import sys
import pymysql

# EnsEMBL connection
# TODO make this prettier
conn = None

def p(line, end=os.linesep):
    """print pretty

    Args:
        line (str): line to print to STDOUT

    Returns:
        pass
    """
    print('\033[93m', '>>> ', line, '\033[0m', end=end)

def query():
    """Queries EnsEMBL. Parameters are HGNC_ID and/or Ensembl_gene_id, whatever is available. Data from EnsEMBLdb will overwrite the client data.
    A(n) identifier(s) should yield one result from EnsEMBLdb. It will be reported if a(n) identifier(s) don't yield any or multiple results.

    Args:
        data (list of dicts): representing the lines and columns in a gene list. The keys of the dicts must match the column names of the EnsEMBLdb query.
        try_hgnc_again (bool): when providing multiple HGNC ids, try until you find a match on EnsEMBLdb. Only one HGNC ID will be used in final result.
    Yields:
        dict: a row with data from ensEMBLdb filled in.
    """
    conn = pymysql.connect(host='ensembldb.ensembl.org', port=5306, user='anonymous', db='homo_sapiens_core_75_37')
    cur = conn.cursor(pymysql.cursors.DictCursor)

    base_query = "select seq_region.name AS Chromosome, g.seq_region_start AS Gene_start, g.seq_region_end AS Gene_stop, g.stable_id AS Ensembl_ID, xg.display_label AS HGNC_ID, g.description, t.stable_id AS Transcript_ID, x.display_label AS RefSeq_ID from gene g left join transcript t on t.gene_id = g.gene_id join xref xg on xg.xref_id = g.display_xref_id join xref x on x.xref_id = t.display_xref_id join seq_region on g.seq_region_id = seq_region.seq_region_id where length(seq_region.name) < 3"

    cur.execute(base_query)
    rs = cur.fetchall() # result set

    Ensembl_ID = ''
    line = ''
    prev_description = ''
    transcripts = []
    print('Chr	Start	Stop	Ensembl_GENE_ID HGNC_Symbol	Ensembl_transcript_to_refseq_transcript	Gene_description')
    for row in rs:
        if row['Ensembl_ID'] != Ensembl_ID:

            if len(transcripts) == 0:
                p(Ensembl_ID + ' has no transcripts!')

            line += '%s:%s' % (Ensembl_ID, '|'.join(transcripts))
            line += '\t' + prev_description if prev_description != None else '';

            print(line)

            # reset
            transcripts = []
            Ensembl_ID = row['Ensembl_ID']

            # sanity check
            Gene_start = row['Gene_start']
            Gene_stop  = row['Gene_stop']
            if row['Gene_start'] > row['Gene_stop']:
                Gene_start, Gene_stop = Gene_stop, Gene_start

            line = '%s\t%d\t%d\t%s\t%s\t' % (row['Chromosome'], Gene_start, Gene_stop, row['Ensembl_ID'], row['HGNC_ID'])
            prev_description = row['description']

        if row['RefSeq_ID'] == None:
            p('%s:%s has no RefSeqID' % (Ensembl_ID, row['Transcript_ID']))

        transcripts.append('%s>%s' % (row['Transcript_ID'], row['RefSeq_ID']))

    # print last one
    line += '|'.join(transcripts)
    line += '\t' + prev_description if prev_description != None else '';
    print(line)

def main(argv):

    # fill in missing blanks
    global conn
    conn = pymysql.connect(host='ensembldb.ensembl.org', port=5306, user='anonymous', db='homo_sapiens_core_75_37')
    data = query()

if __name__ == '__main__':
    main(sys.argv[1:])
