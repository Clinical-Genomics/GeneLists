#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function
import sys
import pymysql
import argparse
import re
import os

gl_header=['Chromosome', 'Gene_start', 'Gene_stop', 'HGNC_ID', 'Disease_group_pathway', 'Protein_name', 'Symptoms', 'Biochemistry', 'Imaging', 'Disease_trivial_name', 'Trivial_name_short', 'Genetic_model', 'OMIM_gene', 'OMIM_morbid', 'Gene_locus', 'Genome_build', 'UniPort_ID', 'Ensembl_gene_id', 'Ensemble_transcript_ID', 'Red_pen', 'Database']

# EnsEMBL connection
# TODO make this prettier
conn = None

# switch: to print or not to print
verbose = False
errors_only = False

def print_header(header=gl_header):
    """Prints the gl_header

    Args:
        header (list, optional): a list of strings
    Returns:
        pass
    Note:
        Prints to STDOUT
    """
    if not errors_only:
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
    if not errors_only:
        ordered_line = list()
        for column_name in gl_header:
            ordered_line.append(str(line[column_name]))
        print("\t".join(ordered_line))

def p(line, end=os.linesep):
    """print only if the verbose switch has been set

    Args:
        line (str): line to print to STDOUT

    Returns:
        pass
    """
    if verbose:
        print(line, end=end)

def query(data, try_hgnc_again=False):
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

    base_query = "select g.seq_region_start AS Gene_start, g.seq_region_end AS Gene_stop, x.display_label AS HGNC_ID, g.stable_id AS Ensembl_gene_id, seq_region.name AS Chromosome from gene g join xref x on x.xref_id = g.display_xref_id join seq_region using (seq_region_id)"
    keys_conds = { 'HGNC_ID': 'x.display_label', 'Ensembl_gene_id': 'g.stable_id', 'Chromosome': 'seq_region.name' }

    keys = ['HGNC_ID', 'Ensembl_gene_id', 'Chromosome'] # these columns will be put into the condition statement if they have a value
    for line in data:
        HGNC_IDs=line['HGNC_ID'].split(',')

        HGNC_ID_i=1
        for HGNC_ID in HGNC_IDs:
            line['HGNC_ID'] = HGNC_ID # actually replace the entry
            conds = [ "%s = %%s" % keys_conds[ key ] for key in keys if key in line and line[key] != None ]
            cond_values = [ line[ key ] for key in keys if key in line ]

            query = "%s where length(seq_region.name) < 3 and %s" % ( base_query, " and ".join(conds) )
            cur.execute(query, cond_values)

            rs = cur.fetchall() # result set

            if len(rs) == 0:
                if HGNC_ID_i == len(HGNC_IDs):
                    not_found_id = HGNC_ID if len(HGNC_IDs) == 1 else HGNC_IDs
                    p("Not found: %s, chromosome: %s" % (not_found_id, line['Chromosome']))
                if not try_hgnc_again: break
            elif len(rs) > 1:
                if HGNC_ID_i > 1:
                    p("Took %s/%s" % (HGNC_ID, HGNC_IDs))
                # this happens when multiple genes are overlapping with a HGNC ID
                p("Multiple entries: %s, chromosome: %s => " % tuple(cond_values), end='')
                p("Adding: %s" % ', '.join(( line['Ensembl_gene_id'] for line in rs )) )
                for entry in rs:
                    yield merge_line(entry, line)
                break
            else:
                if HGNC_ID_i > 1:
                    p("Took %s/%s" % (HGNC_ID, HGNC_IDs))
                for entry in rs:
                    yield merge_line(entry, line)
                break
            HGNC_ID_i += 1

def get_transcript(start, end, ensembl_gene_id=None, hgnc_id=None):
    """Queries EnsEMBL. Parameters are HGNC_ID and/or Ensembl_gene_id, whatever is available. It will return one hit with the ensembl trasncript id.

    Args:
        ensembl_gene_id (str, optional): ensembl_gene_id and/or hgnc_id should be provided.
        hgnc_id (str, optional): ensembl_gene_id and/or hgnc_id should be provided.
        start (int): start coordinate of the possible transcript of this gene.
        end (int): stop coordinate of the possible transcript of this gene.

    Yields:
        dict: with following keys: ensembl_gene_id, hgnc_id, start, end, transcript_id
    """
    cur = conn.cursor(pymysql.cursors.DictCursor)

    base_query = "select g.seq_region_start AS Gene_start, g.seq_region_end AS Gene_stop, x.display_label AS HGNC_ID, g.stable_id AS Ensembl_gene_id, seq_region.name AS Chromosome, t.seq_region_start AS Transcript_start, t.seq_region_end AS Transcript_stop, t.stable_id AS Transcript_id from gene g join xref x on x.xref_id = g.display_xref_id join seq_region using (seq_region_id) join transcript t using (gene_id)"

    conds = {'t.seq_region_start': start, 't.seq_region_end': end}
    if ensembl_gene_id != None:
        conds.update({'g.stable_id': ensembl_gene_id})
    if hgnc_id != None:
        conds.update({'x.display_label': hgnc_id})

    query = "%s where %s" % ( base_query, " and ".join([ '%s = %%s' % column for column in conds.keys() ]) )
    cur.execute(query, [ str(value) for value in conds.values() ])

    rs = cur.fetchall() # result set

    # O-oh .. for now this still means manual intervention!
    if len(rs) > 1:
        p("Getting '%s' ... " % conds.values())
        p('Multiple entries found!')
    elif len(rs) == 0:
        return None
    return rs[0]

def merge_line(ens, client):
    """Will merge dict ens (EnsEMBL data) with client (data). ens will take precedence over client. Changes will be reported.

    Args:
        ens (dict): dict with values for one line of a gene list
        client (dict): with values for one line of a gene list
    Yields:
        dict: merged ens and client dict

    """

    for key, value in client.items():
        if key in ('Gene_start', 'Gene_stop'):
            pass
        elif key in ens:
            if str(ens[key]) != str(value):
                p("%s > %s: ens '%s' diff from client '%s'" % (ens['Ensembl_gene_id'], key, ens[key], value))
        #    else:
        #        p("%s: ens '%s' eq to client '%s'" % (key, ens[key], value))
        #else:
        #    p("%s not in ens!" % key)

    # Check the Gene_start and Gene_stop for being Transcript coordinates
    if 'Gene_start' in client and 'Gene_stop' in client:
        if str(ens['Gene_start']) != str(client['Gene_start']) or str(ens['Gene_stop']) != str(client['Gene_stop']):
            # get the transcript, compare those coordinates and report
            transcript = get_transcript(client['Gene_start'], client['Gene_stop'], ens['Ensembl_gene_id'], ens['HGNC_ID'])
            for key in ('Gene_start', 'Gene_stop'):
                if str(ens[key]) != str(client[key]):
                    if transcript and len(transcript) > 1:
                        p("%s > %s: ens '%s' diff from client '%s', but matches Transcript %s: %s" % (ens['Ensembl_gene_id'], key, ens[key], client[key], transcript['Transcript_id'], transcript[key.replace('Gene', 'Transcript')]))
                    else:
                        p("%s > %s: ens '%s' diff from client '%s'" % (ens['Ensembl_gene_id'], key, ens[key], client[key]))

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

def zero2one(data):
    """Fix 0-based coordinates to 1-based coordinates

    Args:
        data (list of lists): Inner list represents a row in a gene list
    Yields:
        dict: the fixed coordinates data dict
    """
    for line in data:
        for key in ('Gene_start', 'Gene_stop'):
            line[key] = int(line[key]) + 1
        yield line

def download_mim2gene():
    """Download the mim2gene.txt file from omim ftp server. By default the file
    will be downloaded to the location of the script.


    Returns (str): full path and filename of the mim2gene.txt
    """
    from urllib.request import urlretrieve

    filename = os.path.dirname(__file__) + os.path.sep + 'mim2gene.txt'
    (dl_filename, headers) = urlretrieve('ftp://anonymous:kennybilliau%40scilifelab.se@ftp.omim.org/OMIM/mim2gene.txt', filename)
    return dl_filename

def main(argv):
    # set up the argparser
    parser = argparse.ArgumentParser(description='Queries EnsEMBL and fills in the blanks of a gene list. Only columns headers found in gl_headers will be used')
    parser.add_argument('infile', type=argparse.FileType('r'), help='the tsv file with correct headers')
    parser.add_argument('--zero', default=False, action='store_true', dest='zero_based', help='if set, will convert 0-based coordinates to 1-based')
    parser.add_argument('--verbose', default=False, action='store_true', dest='verbose', help='if set, will show conflict messages from EnsEMBLdb inbetween the gene list lines')
    parser.add_argument('--errors-only', default=False, action='store_true', dest='errors_only', help='if set, will not output the gene list, but only the conflict messages from EnsEMBLdb.')
    parser.add_argument('--download-mim2gene', default=False, action='store_true', dest='download_mim2gene', help='if set, will download a new version of the mim2gene.txt file, used to check the OMIM type')
    args = parser.parse_args(argv)

    global verbose
    # make sure we print if we are asked to
    if args.verbose:
        verbose = True

    # show only the EnsEMBLdb conflicts - so verbose, but supress printing of the gene list
    if args.errors_only:
        global errors_only
        verbose = True
        errors_only = True

    # download a new version of mim2gene.txt
    if args.download_mim2gene:
        p('Downloading mim2gene.txt ... ', end='')
        dl_filename = download_mim2gene()
        p('Done: %s' % dl_filename)

    # read in the TSV file
    tsvfile = args.infile
    raw_data = ( line.strip() for line in tsvfile ) # sluuuurp
    parsable_data = ( line.split("\t") for line in raw_data )

    # clean up the input
    clean_data = cleanup(parsable_data)

    # list to dict
    header = next(clean_data) # get the header
    dict_data = list2dict(header, clean_data)

    fixed_data = None
    if args.zero_based:
        # fix 0-based coordinates to be 1-based
        fixed_data = zero2one(dict_data)
    else:
        fixed_data = dict_data

    # fill in missing blanks
    global conn
    conn = pymysql.connect(host='ensembldb.ensembl.org', port=5306, user='anonymous', db='homo_sapiens_core_75_37')
    ensembld_data = query(fixed_data, try_hgnc_again=True)

    # fill in missing values with #NA
    completed_data = fill(ensembld_data)

    # clean up the data from EnsEMBL a bit
    munged_data = munge(completed_data)

    # print the gene list
    print_header()
    for line in munged_data:
        print_line(line)

if __name__ == '__main__':
    main(sys.argv[1:])
