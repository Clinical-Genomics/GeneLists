#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function
import os
import sys
import pymysql
import re

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

def fill_line(row):
    """Formats a line in the regions file with information found in the row

    Args:
        row (dict): with following keys: Chromosome, Gene_start, Gene_stop, Ensembl_ID, HGNC_symbol

    Returns: a tab delimited string

    """

    # sanity check
    Gene_start = row['Gene_start']
    Gene_stop  = row['Gene_stop']
    if row['Gene_start'] > row['Gene_stop']:
        Gene_start, Gene_stop = Gene_stop, Gene_start

    return '%s\t%d\t%d\t%s\t%s\t' % (row['Chromosome'], Gene_start, Gene_stop, row['Ensembl_ID'], row['HGNC_symbol'])

def cleanup_description(description):
    """Remove the comment in the description and clean up invalid characters: ,:;|>

    Args:
        description (str): text to clean up

    Returns: str or None

    """
    if description:
        description = description.strip()
        description = re.sub(r'\[.*\]', '', description)
        description = re.sub(r'[,:;>| ]', '_', description)
        if description.endswith('_'):
            description = description[:-1]
        return description
    return None

def query():
    """Queries EnsEMBL for all transcripts.

    yields (str):
        a line properly BED-formatted
        
    """
    conn = pymysql.connect(host='ensembldb.ensembl.org', port=5306, user='anonymous', db='homo_sapiens_core_75_37')
    cur = conn.cursor(pymysql.cursors.DictCursor)

    """
    external_db_id = 1801
    select * from xref where display_label like 'NM\_%' limit 10;
    """

    base_query = """
        SELECT seq_region.name AS Chromosome, g.seq_region_start AS Gene_start, g.seq_region_end AS Gene_stop, g.stable_id AS Ensembl_ID, xg.display_label AS HGNC_symbol, g.description, t.stable_id AS Transcript_ID, x.dbprimary_acc AS RefSeq_ID
        FROM gene g
        LEFT JOIN transcript t ON t.gene_id = g.gene_id
        JOIN xref xg ON xg.xref_id = g.display_xref_id
        LEFT JOIN object_xref ox ON ox.ensembl_id = t.transcript_id
        JOIN xref x ON x.xref_id = ox.xref_id
        JOIN seq_region ON g.seq_region_id = seq_region.seq_region_id
        WHERE length(seq_region.name) < 3
        AND x.external_db_id in (1801, 1810)
    """

    cur.execute(base_query)
    rs = cur.fetchall() # result set

    row = rs.pop(0)

    # init
    Ensembl_ID = row['Ensembl_ID']
    line = fill_line(row)
    prev_description = cleanup_description(row['description'])
    transcripts = ['%s>%s' % (row['Transcript_ID'], row['RefSeq_ID'])]

    for row in rs:
        if row['Ensembl_ID'] != Ensembl_ID:

            if len(transcripts) == 0:
                p(Ensembl_ID + ' has no transcripts!')

            line += '%s:%s' % (Ensembl_ID, '|'.join(transcripts))
            line += '\t' + prev_description if prev_description != None else '';

            yield line

            # reset
            transcripts = []
            Ensembl_ID = row['Ensembl_ID']

            line = fill_line(row)
            prev_description = cleanup_description(row['description'])

        if row['RefSeq_ID'] == None:
            p('%s:%s has no RefSeqID' % (Ensembl_ID, row['Transcript_ID']))

        transcripts.append('%s>%s' % (row['Transcript_ID'], row['RefSeq_ID']))

    # print last one
    line += '%s:%s' % (Ensembl_ID, '|'.join(transcripts))
    line += '\t' + prev_description if prev_description != None else '';
    yield line

def main(argv):

    # fill in missing blanks
    global conn
    conn = pymysql.connect(host='ensembldb.ensembl.org', port=5306, user='anonymous', db='homo_sapiens_core_75_37')
    data = sorted(query())

    print('#Chromosome	Gene_start	Gene_stop	Ensembl_gene_id	HGNC_symbol	Ensembl_transcript_to_refseq_transcript	Gene_description')
    for line in data:
        print(line)

if __name__ == '__main__':
    main(sys.argv[1:])
